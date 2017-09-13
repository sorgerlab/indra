from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import csv
import time
import tarfile
import zlib
import logging
import xml.etree.ElementTree as ET
import re
import os
import multiprocessing as mp
from os import path
from ftplib import FTP
from io import BytesIO
from indra.util import zip_string, unzip_string
from indra.util import UnicodeXMLTreeBuilder as UTB
from indra.literature.pmc_client import id_lookup
from indra.literature import pubmed_client
from indra.db import get_aws_db, texttypes

logger = logging.getLogger('update_db')

ftp_blocksize = 33554432 # Chunk size recommended by NCBI
BATCH_SIZE = 10000

    


class NihFtpClient(object):
    '''
    High level access to the nih ftp repositories.
    
    Parameters
    ----------
    my_path: str 
        The path to the subdirectory around which this client operates.
    ftp_url: str: (default: 'ftp.ncbi.nlm.nih.gov')
        The url to the ftp site. May be a local directory (see `local`)
    local: bool: (default: False)
        These methods may be run on a local directory (intended for testing).
    '''
    def __init__(self, my_path, ftp_url='ftp.ncbi.nlm.nih.gov', local=False):
        self.my_path = my_path
        self.is_local = local
        self.ftp_url = ftp_url


    def path_join(self, *args):
        joined_str = path.join(*args)
        part_list = joined_str.split('/')
        for part in part_list[1:]:
            if part == '..':
                idx = part_list.index(part) - 1
                part_list.pop(idx)
                part_list.pop(idx)
        return path.join(*part_list)


    def get_ftp_connection(self, ftp_path = None):
        if ftp_path is None:
            ftp_path = self.my_path
        # Get an FTP connection
        ftp = FTP(self.ftp_url)
        ftp.login()
        # Change to the manuscripts directory
        ftp.cwd(ftp_path)
        return ftp


    def get_xml_file(self, xml_file):
        print("Downloading %s" % (xml_file))
        ret = self.get_file(xml_file)
        print("Unzipping")
        xml_bytes = zlib.decompress(ret, 16+zlib.MAX_WBITS)
        print("Parsing XML metadata")
        return ET.XML(xml_bytes, parser=UTB())


    def get_csv_as_dict(self, csv_file):
        csv_str = self.get_file(csv_file).decode('utf8')
        lst = []
        reader = csv.reader(csv_str.splitlines())
        for row in reader:
            lst.append(row)
        return [dict(zip(lst[0], row)) for row in lst[1:]]


    def ret_file(self, f_path, buf):
        full_path = self.path_join(self.my_path, f_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                ftp.retrbinary('RETR /%s' % full_path,
                               callback=lambda s: buf.write(s),
                               blocksize=ftp_blocksize)
                buf.flush()
        else:
            with open(self.path_join(self.ftp_url, full_path), 'rb') as f:
                buf.write(f.read())
                buf.flush()
        return


    def download_file(self, f_path):
        name = path.basename(f_path)
        with open(name, 'wb') as gzf:
            self.ret_file(f_path, gzf)
        return name


    def get_file(self, f_path):
        gzf_bytes = BytesIO()
        self.ret_file(f_path, gzf_bytes)
        return gzf_bytes.getvalue()


    def ftp_ls(self, ftp_path = None):
        if ftp_path is None:
            ftp_path = self.my_path
        else:
            ftp_path = self.path_join(self.my_path, ftp_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                contents = ftp.nlst()
        else:
            contents = os.listdir(self.path_join(self.ftp_url, ftp_path))
        return contents


class Uploader(object):
    '''
    This abstract class provides the api required for any object that is used 
    to upload content into the database.
    '''
    my_source=NotImplemented
    def populate(self, db):
        "A stub for the method used to initially populate the database."
        raise NotImplementedError(
            "`Populate` not implemented for `%s`." % self.__class__.__name__
            )
    
    def update(self, db):
        "A stub for the method used to update the content on the database."
        raise NotImplementedError(
            "`Update` not implemented for `%s`." % self.__class__.__name__
            )

class NihUploader(Uploader):
    '''
    This is an abstract class for all the uploaders that use the nih ftp 
    service.
    
    See `NihFtpClient` for parameters. 
    '''
    my_path = NotImplemented
    def __init__(self, *args, **kwargs):
        self.ftp = NihFtpClient(self.my_path, *args, **kwargs)


class Medline(NihUploader):
    my_path='pubmed/baseline'
    my_source='pubmed'
    def get_deleted_pmids(self):
        del_pmid_str = self.ftp.get_file('../deleted.pmids.gz')
        pmid_list = [
            line.strip() for line in unzip_string(del_pmid_str).split('\n')
            ]
        return pmid_list


    def get_file_list(self):
        all_files = self.ftp.ftp_ls()
        return [k for k in all_files if k.endswith('.xml.gz')]


    def get_article_info(self, xml_file, q = None):
        tree = self.ftp.get_xml_file(xml_file)
        article_info = pubmed_client.get_metadata_from_xml_tree(
            tree, 
            get_abstracts=True, 
            prepend_title=False
            )
        if q is not None:
            q.put(article_info)
            return
        else:
            return article_info

    
    def fix_doi(self, doi):
        "Sometimes the doi is doubled (no idea why). Fix it."
        if doi is None:
            return doi
        L = len(doi)
        if L%2 is not 0:
            return
        if doi[:L] != doi[L:]:
            return
        print("\tFixing doubled doi: %s" % doi)
        return doi[:L]

    def upload_article(self, db, article_info):
        deleted_pmids = self.get_deleted_pmids()
        
        print("%d PMIDs in XML dataset" % len(article_info))
        # Convert the article_info into a list of tuples for insertion into
        # the text_ref table
        text_ref_records = []
        text_content_info = {}
        valid_pmids = set(article_info.keys()).difference(set(deleted_pmids))
        print("%d valid PMIDs" % len(valid_pmids))
        existing_pmids = set(db.get_pmids(valid_pmids))
        print("%d valid PMIDs already in text_refs." % len(existing_pmids))
        pmids_to_add = valid_pmids.difference(existing_pmids)
        print("%d PMIDs to add to text_refs" % len(pmids_to_add))
        for pmid in pmids_to_add:
            pmid_data = article_info[pmid]
            rec = (
                pmid, pmid_data.get('pmcid'), 
                self.fix_doi(pmid_data.get('doi')),
                pmid_data.get('pii')
                )
            text_ref_records.append(
                tuple([None if not r else r for r in rec])
                )
            abstract = pmid_data.get('abstract')
            # Make sure it's not an empty or whitespace-only string
            if abstract and abstract.strip():
                abstract_gz = zip_string(abstract)
                text_content_info[pmid] = \
                    (self.my_source, 'text', texttypes.ABSTRACT, abstract_gz)

        db.copy(
            'text_ref', 
            text_ref_records, 
            ('pmid', 'pmcid', 'doi', 'pii',)
            )

        # Build a dict mapping PMIDs to text_ref IDs
        pmid_list = list(text_content_info.keys())
        tref_list = db.select(
            'text_ref', 
            db.TextRef.pmid.in_([p for p in pmid_list])
            )
        pmid_tr_dict = {pmid:trid for (pmid, trid) in 
                        db.get_values(tref_list, ['pmid', 'id'])}

        # Add the text_ref IDs to the content to be inserted
        text_content_records = []
        for pmid, tc_data in text_content_info.items():
            if pmid not in pmid_tr_dict.keys():
                continue
            tr_id = pmid_tr_dict[pmid] 
            text_content_records.append((tr_id,) + tc_data)


        db.copy(
            'text_content', 
            text_content_records,
            cols=('text_ref_id', 'source', 'format', 'text_type', 'content',)
            )
        return True


    def populate(self, db, n_procs = 1, continuing=False):
        """
        Perform the initial population of the pubmed content into the database.
        
        Inputs
        ------
        db: indra.db.DatabaseManager instance
            The database to which the data will be uploaded.
        n_procs: int
            The number of processes to use when parsing xmls.
        continuing: bool
            If true, assume that we are picking up after an error, or otherwise
            continuing from an earlier process. This means we will skip over
            source files contained in the database. If false, all files will be
            read and parsed.
        """
        xml_files = self.get_file_list()
        tr_list = db.select(
            'source_file', 
            db.SourceFile.source==self.my_source
            )
        existing_files = [tr.name for tr in tr_list]
        
        q = mp.Queue()
        proc_list = []
        for xml_file in xml_files:
            if continuing and xml_file in existing_files:
                print("Skipping %s. Already uploaded." % xml_file)
                continue
            if xml_file not in existing_files:
                db.insert('source_file', source=self.my_source, name=xml_file)
            p = mp.Process(
                target=self.get_article_info,
                args=(xml_file, q, ),
                daemon=True
                )
            proc_list.append(p)
        n_tot = len(proc_list)
        
        for _ in range(n_procs):
            if len(proc_list):
                proc_list.pop(0).start()

        while len(proc_list):
            article_info = q.get() # This will block until at least one is done
            proc_list.pop(0).start()
            self.upload_article(db, article_info)
            n_tot -= 1

        while n_tot is not 0:
            article_info = q.get()
            self.upload_article(db, article_info)
            n_tot -= 1
        
        return


class PmcUploader(NihUploader):
    my_source = NotImplemented
    def __init__(self, *args, **kwargs):
        super(PmcUploader, self).__init__(*args, **kwargs)
        self.tr_cols = ('pmid', 'pmcid', 'doi', 'manuscript_id',)
        self.tc_cols = ('text_ref_id', 'source', 'format', 'text_type', 'content',)


    def upload_batch(self, db, tr_data, tc_data):
        # Get the latest on what is already in the database.
        "Add a batch of text refs and text content to the database."
        pmid_list = [entry['pmid'] for entry in tr_data]
        pmcid_list = [entry['pmcid'] for entry in tr_data]
        #text_refs = db.select(
        #    'text_ref', 
        #    db.TextRef.pmid.in_(pmid_list) | db.TextRef.pmcid.in_(pmcid_list)
        #    )
        #raw_db_conts = db.get_values(
        #    text_refs, 
        #    ['pmid', 'pmcid'],
        #    keyed=True
        #    )
        #del text_refs
        db_conts = {
            'pmid':[tr.pmid for tr in db.select('text_ref', db.TextRef.pmid.in_(pmid_list))],
            'pmcid':[tr.pmcid for tr in db.select('text_ref', db.TextRef.pmcid.in_(pmcid_list))]
            }

        # Define some helpful functions.
        filtered_tr_records = []
        def add_record(tr_entry):
            entry_lst = []
            for k in self.tr_cols:
                if tr_entry[k] is None:
                    entry = tr_entry[k]
                else:
                    entry = tr_entry[k]
                entry_lst.append(entry)
            filtered_tr_records.append(tuple(entry_lst))

        def update_record_with_pmcid(tr_entry):
            tr = db.select_one('text_ref', db.TextRef.pmid==tr_entry['pmid'])
            tr.pmcid = tr_entry['pmcid']
            db.commit('Did not update pmcid %s.' % tr_entry['pmcid'])
        
        def update_record_with_pmid(tr_entry, pmid=None):
            tr = db.select_one('text_ref', db.TextRef.pmcid==tr_entry['pmcid'])
            tr.pmid = tr_entry['pmid'] if pmid is None else pmid
            db.commit('Failed to update pmid %s.' % tr.pmid)

        # Process the text ref data.
        for tr_entry in tr_data:
            # Try to check by pmid first
            if tr_entry['pmid'] is None:
                pmid = id_lookup(tr_entry['pmcid'])['pmid']
                if pmid is None:
                    if tr_entry['pmcid'] not in db_conts['pmcid']:
                        # Add a record without a pmid (it happens)
                        add_record(tr_entry)
                elif pmid not in db_conts['pmid']:
                    if tr_entry['pmcid'] not in db_conts['pmcid']:
                        tr_entry['pmid'] = pmid
                        add_record(tr_entry)
                    else:
                        update_record_with_pmid(tr_entry, pmid)
                else:
                    if tr_entry['pmcid'] not in db_conts['pmcid']:
                        update_record_with_pmcid(tr_entry)
            elif tr_entry['pmid'] not in db_conts['pmid']:
                if tr_entry['pmcid'] not in db_conts['pmcid']:
                    add_record(tr_entry)
                else:
                    update_record_with_pmid(tr_entry)
            else:
                if tr_entry['pmcid'] not in db_conts['pmcid']:
                    update_record_with_pmcid(tr_entry)
        
        # Upload the text content data.
        print('Adding %d new text refs...' % len(filtered_tr_records))
        db.copy('text_ref', filtered_tr_records, self.tr_cols)

        # Process the text content data
        arc_pmcid_list = [
            tr['pmcid' ]
            for tr in tr_data
            ]
        tref_list = db.select(
            'text_ref', 
            db.TextRef.pmcid.in_(arc_pmcid_list)
            )
        pmcid_trid_dict = {
            pmcid:trid for (pmcid, trid) in
            db.get_values(tref_list, ['pmcid', 'id'])
            }
        existing_tcs = db.select(
            'text_content', 
            db.TextContent.text_ref_id.in_(pmcid_trid_dict.values()),
            db.TextContent.source==self.my_source,
            db.TextContent.format=='xml'
            ) # This should be a very small list, in general.
        existing_tc_records = [
            (tc.text_ref_id, tc.source, tc.format, tc.text_type)
            for tc in existing_tcs
            ]
        tc_records = [
            (pmcid_trid_dict[pmcid], self.my_source, 'xml', txt_type, cont)
            for pmcid, txt_type, cont in tc_data
            ]
        filtered_tc_records = [
            rec for rec in tc_records if rec[:-1] not in existing_tc_records
            ]

        # Upload the text content data.
        print('Adding %d more text content entries...' % len(filtered_tc_records))
        db.copy('text_content', filtered_tc_records, self.tc_cols)


    def upload_archive(self, db, archive, q=None):
        "Process a single tar gzipped article archive."
        tr_data = []
        tc_data = []
        with tarfile.open(archive, mode='r:gz') as tar:
            print('Loading...')
            xml_files = [m for m in tar.getmembers() if m.isfile()]
            for i, xml_file in enumerate(xml_files):
                #print("Reading %s. File %d/%d." % (xml_file, i, len(xml_files)))
                xml_str = tar.extractfile(xml_file).read()
                try:
                    tree = ET.XML(xml_str)
                except ET.ParseError:
                    print("Could not parse %s. Skipping." % xml_file.name)
                    continue
                id_data = {
                    e.get('pub-id-type'):e.text for e in 
                    tree.findall('.//article-id')
                    }
                if 'pmc' not in id_data.keys():
                    print("Did not get a 'pmc' in %s." % xml_file)
                    continue
                    
                #if 'PMC'+id_data['pmc'] not in pmcid_list:
                #    continue
                if 'pmcid' not in id_data.keys():
                    id_data['pmcid'] = 'PMC' + id_data['pmc']
                rec_dict = dict.fromkeys(self.tr_cols)
                rec_dict.update(id_data)
                tr_data.append(rec_dict)
                for typ, lbl in [(texttypes.ABSTRACT,'abstract'), (texttypes.FULLTEXT,'body')]:
                    cont_xml = tree.find('.//'+lbl)
                    if cont_xml is not None:
                        content = ET.tostring(cont_xml).decode('utf8')
                        tc_data.append(
                            (
                                id_data['pmcid'], 
                                typ, 
                                zip_string(content)
                                )
                            )
                if (i+1)%BATCH_SIZE==0: # Upload in batches, so as not to overwhelm ram.
                    print("Submitting batch %d of data for %s..." % 
                          ((i+1)/BATCH_SIZE, archive))
                    if q is not None:
                        q.put((tr_data, tc_data))
                    else:
                        self.upload_batch(db, tr_data, tc_data)
                    tr_data = []
                    tc_data = []
            else:
                print("Submitting final batch of data for %s..." % archive)
                if q is not None:
                    q.put((tr_data, tc_data))
                else:
                    self.upload_batch(db, tr_data, tc_data)
        return


    def process_archive(self, db, archive, q=None):
        try:
            print('Downloading archive %s.' % archive)
            self.ftp.download_file(archive)
            self.upload_archive(db, archive, q)
        finally:
            os.remove(archive)


    def get_file_list(self):
        return [k for k in self.ftp.ftp_ls() if self.is_archive(k)]
    
    
    def populate(self, db, n_procs=1, continuing=False):
        """
        Perform the initial population of the pubmed content into the database.
        
        Inputs
        ------
        db: indra.db.DatabaseManager instance
            The database to which the data will be uploaded.
        n_procs: int
            (Not used) The number of processes to use when parsing xmls.
        continuing: bool
            If true, assume that we are picking up after an error, or
            otherwise continuing from an earlier process. This means we will
            skip over source files contained in the database. If false, all
            files will be read and parsed.
        """
        archives = self.get_file_list()
        
        tr_list = db.select(
            'source_file', 
            db.SourceFile.source==self.my_source
            )
        existing_arcs = [tr.name for tr in tr_list]
        
        q = mp.Queue(len(archives))
        wait_list = []
        for archive in archives:
            if continuing and archive in existing_arcs:
                print("Skipping %s. Already uploaded." % archive)
                continue
            p = mp.Process(
                target=self.process_archive,
                args=(db,archive,q,),
                daemon=True
                )
            wait_list.append((archive, p))
        
        active_list = []
        def start_next_proc():
            if len(wait_list) is not 0:
                archive, proc = wait_list.pop(0)
                proc.start()
                active_list.append((archive, proc))
        
        # Start the processes running
        for _ in range(min(n_procs, len(archives))):
            start_next_proc()

        # Monitor the processes while any are still active.
        while len(active_list) is not 0:
            for a, p in [(a,p) for a,p in active_list if not p.is_alive()]:
                if a not in existing_arcs:
                    db.insert('source_file', source=self.my_source, name=a)
                active_list.remove((a, p))
                start_next_proc()
            try:
                tr_data, tc_data = q.get_nowait() # This will block until at least one is done
            except:
                continue
            self.upload_batch(db, tr_data, tc_data)
            time.sleep(0.1)
                
        # Empty the queue.
        while not q.empty():
            try:
                tr_data, tc_data = q.get(timeout=1)
            except:
                break
            self.upload_batch(db, tr_data, tc_data)
        
        return


class PmcOA(PmcUploader):
    my_path = 'pub/pmc'
    my_source = 'pmc_oa'
    
    def is_archive(self, k):
        return k.startswith('articles') and k.endswith('.xml.tar.gz')


class Manuscripts(PmcUploader):
    my_path = 'pub/pmc/manuscript'
    my_source = 'manuscripts'
    def get_tarname_from_filename(self, fname):
        "Get the name of the tar file based on the file name (or a pmcid)."
        re_match = re.match('(PMC00\d).*?', fname)
        if re_match is not None:
            tarname = re_match.group(0) + 6*'X' + '.xml.tar.gz'
        else:
            tarname = None
        return tarname


    def is_archive(self, k):
        return k.endswith('.xml.tar.gz')


if __name__ == '__main__':
    from sys import argv
    if '--continue' in argv or 'continue' in argv:
        continuing=True
    else:
        continuing=False
    if '-n' in argv:
        n = int(argv[argv.index('-n') + 1])
    db = get_aws_db()
    if not continuing:
        db._clear()
    Medline().populate(db, n, continuing)
    PmcOA().populate(db, n, continuing)
    Manuscripts().populate(db, n, continuing)

    # High-level content update procedure
    # 1. Download MEDLINE baseline, will contain all PMIDs, abstracts,
    #    other info
    # 2. Download PMC, update text refs with PMCIDs where possible, with
    #    update commands.
    # 3. Add PMC-OA content.
    # 4. Add PMC-author manuscript information, updating files with manuscript
    #    IDs and content.
    # 5. (Optional): Run script to check for/obtain DOIs for those that are
    #    missing
    # 6. Obtain Elsevier XML for Elsevier articles
    # 7. Obtain Springer content for Springer articles
    # 8. Obtain scraped PDF content for other articles available via OA-DOI.


