from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import csv
import time
import shutil
import tarfile
import tempfile
import functools
import multiprocessing as mp
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

logger = logging.getLogger('update_db')

PmcAuthInfo = namedtuple('PmcAuthInfo', ('File', 'PMCID', 'PMID', 'MID'))
PmcOaInfo = namedtuple('PmcOaInfo',
                       ('File', 'Article_Citation', 'PMCID',
                        'Last_Updated', 'PMID', 'License'))

ftp_blocksize = 33554432 # Chunk size recommended by NCBI


class Progenetor(object):
    ftp_path = NotImplemented
    def __init__(self, db):
        self.db = db

    def get_ftp_connection(self, ftp_path = None):
        if ftp_path is None:
            ftp_path = self.ftp_path
        pmc_ftp_url = 'ftp.ncbi.nlm.nih.gov'
        # Get an FTP connection
        ftp = FTP(pmc_ftp_url)
        ftp.login()
        # Change to the manuscripts directory
        ftp.cwd(ftp_path)
        return ftp


    def get_xml_data(self, xml_file):
        ftp = self.get_ftp_connection()
        #for ix, xml_file in enumerate(xml_files):
        # Download the Gzipped content into a BytesIO
        gzf_bytes = BytesIO()
        #print("Downloading %s (%d of %d)" % (xml_file, ix+1, len(xml_files)))
        print("Downloading %s" % (xml_file))
        ftp.retrbinary('RETR %s' % xml_file,
                                callback=lambda b: gzf_bytes.write(b),
                                blocksize=ftp_blocksize)
        # Unzip the BytesIO into uncompressed bytes
        print("Unzipping")
        xml_bytes = zlib.decompress(gzf_bytes.getvalue(), 16+zlib.MAX_WBITS)
        # Convert the bytes into an XML ElementTree
        print("Parsing XML metadata")
        tree = ET.XML(xml_bytes, parser=UTB())
        ftp.close()
        return tree


    def get_records_and_content(self, xml_file):
        raise NotImplementedError()


    def upload_xml_file(self, xml_file):
        tr_records, tc_records = self.get_records_and_content(xml_file)
        return True
    
    
    def get_xml_file_list(self):
        raise NotImplementedError()
    

    def ftp_ls(self, ftp_path = None):
        if ftp_path is None:
            ftp_path = self.ftp_path
        else:
            ftp_path = path.join(self.ftp_path, ftp_path)
        ftp = self.get_ftp_connection()
        contents = dict(ftp.mlsd())
        ftp.close()
        return contents
    
    
    def populate(self):
        xml_files = self.get_xml_list()
        self.upload_xml_file(xml_files[0])
    

class Medline(Progenetor):
    ftp_path='/pubmed/baseline'
    def get_deleted_pmids(self):
        ftp = self.get_ftp_connection('/pubmed')
        file_gz = BytesIO()
        ftp.retrbinary('RETR deleted.pmids.gz', callback=lambda b: file_gz.write(b),
                       blocksize=ftp_blocksize)
        pmid_list = [line.strip()
                     for line in unzip_string(file_gz.getvalue()).split('\n')]
        ftp.close()
        return pmid_list


    def get_xml_list(self):
        all_files = self.ftp_ls()
        return [k for k in all_files.keys() if k.endswith('.xml.gz')]

    
    def get_records_and_content(self, xml_file):
        deleted_pmids = self.get_deleted_pmids()
        
        tree = self.get_xml_data(xml_file)
        
        # Get the article metadata from the tree
        article_info = pubmed_client.get_metadata_from_xml_tree(
                        tree, get_abstracts=True, prepend_title=False)
        print("%d PMIDs in XML dataset" % len(article_info))
        # Convert the article_info into a list of tuples for insertion into
        # the text_ref table
        text_ref_records = []
        text_content_info = {}
        valid_pmids = set(article_info.keys()).difference(set(deleted_pmids))
        print("%d valid PMIDs" % len(valid_pmids))
        existing_pmids = set(self.db.get_all_pmids())
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
        
        
        self.db.copy(
            'text_ref', 
            text_ref_records, 
            ('pmid', 'pmcid', 'doi', 'pii',)
            )
        
        # Build a dict mapping PMIDs to text_ref IDs
        pmid_list = list(text_content_info.keys())
        tref_list = self.db.select('text_ref', self.db.TextRef.pmid.in_([p.encode() for p in pmid_list]))
        pmid_tr_dict = {pmid.decode('utf-8'):trid for (pmid, trid) in 
                        self.db.get_values(tref_list, ['pmid', 'id'])}
    
        # Add the text_ref IDs to the content to be inserted
        text_content_records = []
        missing_pmids = 0
        for pmid, tc_data in text_content_info.items():
            if pmid not in pmid_tr_dict.keys():
                #print('pmid %s not found.' % pmid)
                missing_pmids += 1
                continue
            tr_id = pmid_tr_dict[pmid] 
            text_content_records.append((tr_id,) + tc_data)
        print("We missed %d/%d pmids." % (missing_pmids, len(text_content_info)))
        
        
        self.db.copy(
            'text_content', 
            text_content_records,
            cols=('text_ref_id', 'source', 'format', 'text_type', 'content',)
            )
        
        return text_ref_records, text_content_records


def initialize_medline(db):
    Medline(db).populate()

class PmcOA(Progenetor):
    ftp_path = '/pub/pmc'
    def get_file_info(self, filename, datatype):
        ftp = self.get_ftp_connection(self.ftp_path)
        # Get the list of files from the tab-separated .txt file
        print("Downloading %s" % filename)
        filelist_bytes = []
        ftp.retrbinary('RETR %s' % filename,
                      callback=lambda b: filelist_bytes.append(b),
                      blocksize=ftp_blocksize)
        ftp.close()
        # Process the file info (skip the header line)
        print("Processing %s" % filename)
        file_str = StringIO(b''.join(filelist_bytes).decode('ascii'))
        csv_reader = csv.reader(file_str, delimiter=',', quotechar='"')
        next(csv_reader) # Skip the header row
        info_list = [datatype(*row) for row in csv_reader if row]
        return info_list

    def populate(self):
        article_archs = [k for k in self.ftp_ls().keys() if k.startswith('articles') and k.endswith('.xml.tar.gz')]
        current_id_list = self.db.get_values(self.db.select('text_ref'), ['pmid', 'pmcid'])
        pmcid_list = [i[1].decode('utf8') for i in current_id_list if i[1] is not None]
        
        tr_records = []
        tc_records = []
        for article_arch in ['articles.A-B.xml.tar.gz']:#article_archs[:1]:
            try:
                #with GzipFile(article_arch, 'w') as gzf:
                #    ftp = self.get_ftp_connection()
                #    ftp.retrbinary('RETR %s' % article_arch,
                #                   callback=lambda s: gzf.write(s),
                #                   blocksize=ftp_blocksize)
                #    gzf.flush()
                #    ftp.close()
                with GzipFile(article_arch, 'r') as gzf:
                    tar = tarfile.open(fileobj=gzf, mode='r:gz')
                    xml_files = [m for m in tar.getmembers() if m.isfile()]
                    for i, xml_file in enumerate(xml_files):
                        #print("Reading %s. File %d/%d." % (xml_file, i, len(xml_files)))
                        xml_str = tar.extractfile(xml_file).read()
                        try:
                            tree = ET.XML(xml_str)
                        except ET.ParseError:
                            print("Could not parse %s. Skipping." % xml_file)
                            continue
                        id_data = {e.get('pub-id-type'):e.text for e in tree.findall('.//article-id')}
                        abst_xml = tree.find('.//abstract')
                        fulltext_xml = tree.find('.//body')
                        if 'PMC'+id_data['pmc'] not in pmcid_list:
                            continue
                        print("Gonna add %s to db." % 'PMC'+id_data['pmc'])
                        tr_records.append(id_data)
                        tc_records.append((abst_xml, fulltext_xml))
                        if (i+1)%10000==0:
                            print(i)
                            tr_records = []
                            tc_records = []
            finally:
                pass
                #os.remove(article_arch)
        return (tr_records, tc_records)




        


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


def initialize_pmc_oa(db):
    #ftp_path = '/pub/pmc'
    PmcOA(db).populate()
    # The high-level procedure:
    # 1. Get info on all Open Access Subset papers currently in PMC
    #pmc_info_list = _get_file_info(ftp_path, 'oa_file_list.csv', PmcOaInfo)
    # 2. Add text_refs to database for any we don't currently have indexed
    #_update_text_refs(pmc_info_list, source='pmc_oa')
    # 3. Find out which ones are missing pmc_oa_xml content in the databse
    #pmc_info_missing = get_missing_oa_xml_pmcids(pmc_info_list)
    # 4. Figure out which archives we'll need to download
    #xml_list = get_xml_archives_to_download(pmc_info_missing)
    # 5. To save space, download and extract files one at a time
    #xml_list = ['PMC004XXXXXX.xml.tar.gz']
    #for xml_file in xml_list:
    #    tmp_dir = tempfile.mkdtemp(prefix='tmpIndra', dir='.')
    #    download_xml_archive(xml_file, tmp_dir)
    #    update_text_content(pmc_info_list, xml_file, tmp_dir)
    #    #shutil.rmtree(tmp_dir)


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
    download_baseline(db) # Eventually this should only be done once/year.
    #initialize_pmc_oa(db)
    initialize_pmc_manuscripts(db)

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


