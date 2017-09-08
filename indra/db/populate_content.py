from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import csv
import time
import tarfile
import zlib
import logging
import xml.etree.ElementTree as ET
from indra.util import write_unicode_csv, zip_string, unzip_string
from indra.util import UnicodeXMLTreeBuilder as UTB
from ftplib import FTP
from collections import namedtuple
from indra.db import get_aws_db
from io import StringIO, BytesIO
from indra.literature import pubmed_client
from gzip import GzipFile
from os import path
from indra.literature.pmc_client import id_lookup

logger = logging.getLogger('update_db')

PmcAuthInfo = namedtuple('PmcAuthInfo', ('File', 'PMCID', 'PMID', 'MID'))
PmcOaInfo = namedtuple('PmcOaInfo',
                       ('File', 'Article_Citation', 'PMCID',
                        'Last_Updated', 'PMID', 'License'))

ftp_blocksize = 33554432 # Chunk size recommended by NCBI

def path_join(*args):
    joined_str = path.join(*args)
    part_list = joined_str.split('/')
    for part in part_list[1:]:
        if part == '..':
            idx = part_list.index(part) - 1
            part_list.pop(idx)
            part_list.pop(idx)
    return path.join(*part_list)


class Progenetor(object):
    my_path = NotImplemented
    def __init__(self, ftp_url='ftp.ncbi.nlm.nih.gov', local=False):
        self.is_local = local
        self.ftp_url = ftp_url


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
        full_path = path_join(self.my_path, f_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                ftp.retrbinary('RETR /%s' % full_path,
                               callback=lambda s: buf.write(s),
                               blocksize=ftp_blocksize)
                buf.flush()
        else:
            with open(path_join(self.ftp_url, full_path), 'rb') as f:
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
            ftp_path = path_join(self.my_path, ftp_path)
        if not self.is_local:
            with self.get_ftp_connection() as ftp:
                contents = ftp.nlst()
        else:
            contents = os.listdir(path_join(self.ftp_url, ftp_path))
        return contents


    def populate(self, db):
        raise NotImplementedError()


class Medline(Progenetor):
    my_path='pubmed/baseline'
    def get_deleted_pmids(self):
        del_pmid_str = self.get_file('../deleted.pmids.gz')
        pmid_list = [
            line.strip() for line in unzip_string(del_pmid_str).split('\n')
            ]
        return pmid_list


    def get_xml_list(self):
        all_files = self.ftp_ls()
        return [k for k in all_files if k.endswith('.xml.gz')]


    def upload_xml_file(self, db, xml_file):
        deleted_pmids = self.get_deleted_pmids()
        
        tree = self.get_xml_file(xml_file)
        
        # Get the article metadata from the tree
        try:
            article_info = pubmed_client.get_metadata_from_xml_tree(
                        tree, get_abstracts=True, prepend_title=False)
        except Exception as e:
            raise e
        print("%d PMIDs in XML dataset" % len(article_info))
        # Convert the article_info into a list of tuples for insertion into
        # the text_ref table
        text_ref_records = []
        text_content_info = {}
        valid_pmids = set(article_info.keys()).difference(set(deleted_pmids))
        print("%d valid PMIDs" % len(valid_pmids))
        existing_pmids = set(db.get_all_pmids())
        print("%d existing PMIDs in text_refs" % len(existing_pmids))
        pmids_to_add = valid_pmids.difference(existing_pmids)
        print("%d PMIDs to add to text_refs" % len(pmids_to_add))
        for pmid in pmids_to_add:
            pmid_data = article_info[pmid]
            rec = (pmid, pmid_data.get('pmcid'), pmid_data.get('doi'),
                   pmid_data.get('pii'))
            text_ref_records.append(tuple([None if not r else r.encode('utf8')
                                          for r in rec]))
            abstract = pmid_data.get('abstract')
            # Make sure it's not an empty or whitespace-only string
            if abstract and abstract.strip():
                abstract_gz = zip_string(abstract)
                text_content_info[pmid] = \
                                (b'pubmed', b'text', b'abstract', abstract_gz)

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


    def populate(self, db):
        xml_files = self.get_xml_list()
        for xml_file in xml_files:
            self.upload_xml_file(db, xml_file)


class PmcOA(Progenetor):
    my_path = 'pub/pmc'
    def upload_batch(self, db, tr_data, tc_data):
        # Get the latest on what is already in the database.
        "Add a batch of text refs and text content to the database."
        current_id_list = db.get_values(
            db.select('text_ref'), 
            ['pmid', 'pmcid']
            )
        db_conts = {}
        for i, id_type in enumerate(['pmid', 'pmcid']):
            db_conts[id_type] = [
                e[i] for e in current_id_list if e[i] is not None
                ]

        # Define some helpful functions.
        filtered_tr_records = []
        def add_record(tr_entry):
            entry_lst = []
            for k in self.tr_cols:
                if tr_entry[k] is None:
                    entry = tr_entry[k]
                else:
                    entry = tr_entry[k].encode('utf8')
                entry_lst.append(entry)
            filtered_tr_records.append(tuple(entry_lst))

        def update_record(tr_entry):
            tr = db.select('text_ref', db.TextRef.pmid==tr_entry['pmid'])
            tr.pmcid = tr_entry['pmcid']
            db.commit('Did not update pmcid %s.' % tr_entry['pmcid'])

        # Process the text ref data.
        for tr_entry in tr_data:
            # Try to check by pmid first
            if tr_entry['pmid'] is None:
                print("Did not get pmid for %s." % tr_entry['pmcid'])
                pmid = id_lookup(tr_entry['pmcid'])['pmid']
                if pmid is None:
                    if tr_entry['pmcid'] not in db_conts['pmcid']:
                        # Add a record without a pmid (it happens)
                        add_record(tr_entry)
                elif pmid not in db_conts['pmid']:
                    add_record(tr_entry)
                else:
                    if tr_entry['pmcid'] not in db_conts['pmcid']:
                        update_record(tr_entry)
            elif tr_entry['pmid'] not in db_conts['pmid']:
                # If there is not pmid, there is no pmcid in db.
                add_record(tr_entry)
            else:
                if tr_entry['pmcid'] not in db_conts['pmcid']:
                    update_record(tr_entry)

        # Upload the text content data.
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
        pmcid_tr_dict = {
            pmcid:trid for (pmcid, trid) in
            db.get_values(tref_list, ['pmcid', 'id'])
            }
        tc_records = [
            (pmcid_tr_dict[pmcid], 'pmc-oa', 'xml', txt_type, cont)
            for pmcid, txt_type, cont in tc_data
            ]

        # Upload the text content data.
        db.copy('text_content', tc_records, self.tc_cols)


    def upload_article_archive(self, db, article_arch):
        "Process a single tar gzipped article archive."
        tr_data = []
        tc_data = []
        
        self.download_file(article_arch)
        
        with tarfile.open(article_arch, mode='r:gz') as tar:
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
                for typ, lbl in [('abstract','abstract'), ('fulltext','body')]:
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
                if (i+1)%10000==0: # Upload in batches, so as not to overwhelm ram.
                    self.upload_batch(db, tr_data, tc_data)
                    tr_data = []
                    tc_data = []
            else:
                self.upload_batch(db, tr_data, tc_data)
        return


    def populate(self, db):
        self.tr_cols = ('pmid', 'pmcid', 'doi', 'manuscript_id',)
        self.tc_cols = ('text_ref_id', 'source', 'format', 'text_type', 'content',)
        
        article_archs = [
            k for k in self.ftp_ls() 
            if k.startswith('articles') and k.endswith('.xml.tar.gz')
            ]
        for article_arch in article_archs:
            try:
                self.upload_article_archive(db, article_arch)
            except Exception as e:
                raise e
            finally:
                os.remove(article_arch)
        return


def _update_text_refs(pmc_info_list):
    # Insert any missing text_refs into database.
    stored_pmcids = [tr.pmcid for tr in db.select('text_ref', db.TextRef.pmcid!=None)]
    missing_set = set([p.PMCID for p in pmc_info_list]).difference(
                                                      set(stored_pmcids))
    pmc_info_missing = [p for p in pmc_info_list if p.PMCID in missing_set]
    # Write the missing records to a CSV so we can use the copy command
    pmc_data = []
    for pi in pmc_info_missing:
        if pi.PMID == '0' or pi.PMID == '':
            pmid = None
        else:
            pmid = pi.PMID.encode('ascii')
        if isinstance(pi, PmcAuthInfo):
            manuscript_id = pi.MID.encode('ascii')
        else:
            manuscript_id = None
        pmc_data.append((pmid, pi.PMCID.encode('ascii'), manuscript_id))
    # Copy using pgcopy
    db.copy('text_ref', pmc_data, cols=('pmid', 'pmcid', 'manuscript_id'))



class Manuscripts(Progenetor):
    my_path = 'pub/pmc/manuscript'
    def populate(self, db):
        pass

def initialize_pmc_manuscripts(db):
    ftp_path = '/pub/pmc/manuscript'

    def get_missing_auth_xml_pmcids(pmc_info_list):
        # Get list of PMCIDs for which we've already stored author manuscripts
        stored_pmc_ids = db.get_auth_xml_pmcids()
        print("%d pmc_auth_xml PMCIDs already in DB" % len(stored_pmc_ids))
        missing_set = set([p.PMCID for p in pmc_info_list]).difference(
                                                          set(stored_pmc_ids))
        pmc_info_missing = [p for p in pmc_info_list if p.PMCID in missing_set]
        print("%d pmc_auth_xml PMCIDs left to load in DB" %
              len(pmc_info_missing))
        return pmc_info_missing

    def get_xml_archives_to_download(pmc_info_list):
        ftp = _get_ftp_connection(ftp_path)
        # Get the list of .xml.tar.gz files
        xml_files = [f[0] for f in ftp.mlsd() if f[0].endswith('.xml.tar.gz')]
        print("xml_files: %s" % xml_files)
        # Get list of stems in the list of missing files
        stems = set([])
        for pi in pmc_info_list:
            stems.add(pi.File[0:6])
        to_download = [f for f in xml_files if f[0:6] in stems]
        ftp.close()
        return to_download

    def download_xml_archive(filename, tmp_dir):
        ftp = _get_ftp_connection(ftp_path)
        # Function to write to local file with progress updates
        def write_to_file(fp, b):
            fp.write(b)
        outfilepath = os.path.join(tmp_dir, filename)
        print("Getting %s" % filename)
        with open(outfilepath, 'wb') as f:
            ftp.retrbinary('RETR %s' % filename, callback=lambda b: f.write(b),
                           blocksize=ftp_blocksize)
        # Extract all files in the TAR archive
        tf = tarfile.open(outfilepath)
        print("Extracting all files from %s" % outfilepath)
        tf.extractall(path=tmp_dir)
        print("Done extracting files.")
        ftp.close()

    def update_text_content(pmc_info_list, xml_file, tmp_dir, source, poolsize=4):
        # Filter PMCIDs to store to only those contained in this XML file
        stem = xml_file[0:6]
        pmc_to_store = [p for p in pmc_info_list if p.File.startswith(stem)]
        import random
        random.shuffle(pmc_to_store)
        # Get the text ref IDs by PMCID
        pmcid_tr_dict = dict(db.select(
            'text_ref',
            db.TextRef.pmcid.in_([pi.PMCID for pi in pmc_to_store])
            ))
        pmc_blocksize = 2000
        start_ix_list = range(0, len(pmc_to_store), pmc_blocksize)
        for start_ix in start_ix_list:
            content_block_rows = []
            t0 = time.time()
            end_ix = start_ix + pmc_blocksize
            if end_ix > len(pmc_to_store):
                end_ix = len(pmc_to_store)
            for pi in pmc_to_store[start_ix:end_ix]:
                xml_path = os.path.join(tmp_dir, pi.File)
                text_ref_id = pmcid_tr_dict[pi.PMCID].id
                if os.path.exists(xml_path):
                    # Look up the text_ref_id for this PMCID
                    # Read the XML file in text mode
                    with open(xml_path, 'rt') as f:
                        content = f.read()
                    # Compress the content
                    content_gz = zip_string(content)
                    # Add to our CSV rows
                    content_block_rows.append((text_ref_id, source, 'pmc_auth_xml', 'fulltext',
                                               content_gz))
                else:
                    print("Could not find file %s" % xml_path)
            t1 = time.time()
            db.copy('text_content', content_block_rows, ('text_ref_id', 'source', 'format', 'text_type', 'content'))

            print("Copied files %d to %d from %s" %
                  (start_ix, end_ix, xml_file))
            end = time.time()
            read_time = (1000 * (t1 - t0)) / pmc_blocksize
            write_db_time = (1000 * (end - t1)) / pmc_blocksize
            total_time = (1000 * (end - t0)) / pmc_blocksize
            print("Total time: %.1f (read %.1f, write DB %.1f)"
                  % (total_time, read_time, write_db_time))

    # The high-level procedure:
    # 1. Get info on all author manuscripts currently in PMC
    pmc_info_list = _get_file_info(ftp_path, 'filelist.csv', PmcAuthInfo)
    # 2. Add text_refs to database for any we don't currently have indexed
    _update_text_refs(pmc_info_list)
    # 3. Find out which ones are missing auth_xml content in the databse
    pmc_info_missing = get_missing_auth_xml_pmcids(pmc_info_list)
    # 4. Figure out which archives we'll need to download
    #xml_list = get_xml_archives_to_download(pmc_info_missing)
    # 5. To save space, download and extract files one at a time
    xml_list = ['PMC004XXXXXX.xml.tar.gz']
    for xml_file in xml_list:
        #tmp_dir = tempfile.mkdtemp(prefix='tmpIndra', dir='.')
        tmp_dir = 'tmpIndrac2km4wir'
        #download_xml_archive(xml_file, tmp_dir)
        update_text_content(pmc_info_list, xml_file, tmp_dir, 'pmc_auth')
        #shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    db = get_aws_db()
    db._clear()
    Medline().populate(db)
    PmcOA().populate(db)
    Manuscripts().populate(db)

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


