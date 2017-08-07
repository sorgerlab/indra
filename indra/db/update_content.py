import os
from io import StringIO, BytesIO
import csv
import time
import shutil
import tarfile
import tempfile
from ftplib import FTP
from collections import namedtuple
from indra import db
import pgcopy
import functools
import multiprocessing as mp
from indra.util import write_unicode_csv, zip_string


PmcInfo = namedtuple('PmcInfo', ('File', 'PMCID', 'PMID', 'MID'))

ftp_blocksize = 33554432 # Chunk size recommended by NCBI

def initialize_pmc_manuscripts():

    def get_ftp_connection():
        pmc_ftp_url = 'ftp.ncbi.nlm.nih.gov'
        # Get an FTP connection
        ftp = FTP(pmc_ftp_url)
        ftp.login()
        # Change to the manuscripts directory
        ftp.cwd('/pub/pmc/manuscript')
        return ftp

    def get_file_info():
        ftp = get_ftp_connection()
        # Get the list of files from the CSV file
        tmp_dir = tempfile.mkdtemp(prefix='tmpIndra', dir='.')
        filelist_bytes = []
        print("Downloading filelist.csv")
        ftp.retrbinary('RETR filelist.csv',
                      callback=lambda b: filelist_bytes.append(b),
                      blocksize=ftp_blocksize)
        filelist_csv = b''.join(filelist_bytes).decode('ascii').split('\n')
        # Namedtuple for working with PMC info entries
        print("Processing filelist.csv")
        # Process the file info (skip the header line)
        pmc_info_list = [PmcInfo(*line.split(',')) for line in filelist_csv
                         if line][1:]
        shutil.rmtree(tmp_dir)
        ftp.close()
        return pmc_info_list

    def update_text_refs(pmc_info_list):
        # Insert any missing text_refs into database.
        stored_pmc_ids = [r[0] for r in db.select('text_ref', 'pmcid')]
        missing_set = set([p.PMCID for p in pmc_info_list]).difference(
                                                          set(stored_pmc_ids))
        pmc_info_missing = [p for p in pmc_info_list if p.PMCID in missing_set]
        # Write the missing records to a CSV so we can use the copy command
        pmc_info_to_copy = []
        for pi in pmc_info_missing:
            if pi.PMID == '0':
                pmid = None
            else:
                pmid = pi.PMID.encode('ascii')
            pmc_info_to_copy.append((b'pmc', pmid, pi.PMCID.encode('ascii'),
                                     pi.MID.encode('ascii')))
        # Copy using pgcopy
        conn = db.get_connection()
        cols = ('source', 'pmid', 'pmcid', 'manuscript_id')
        mgr = pgcopy.CopyManager(conn, 'text_ref', cols)
        mgr.copy(pmc_info_to_copy, BytesIO)

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
        ftp = get_ftp_connection()
        # Get the list of .xml.tar.gz files
        xml_files = [f[0] for f in ftp.mlsd() if f[0].endswith('.xml.tar.gz')]
        print("xml_files: %s" % xml_files)
        # Get list of stems in the list of missing files
        stems = set([])
        for pi in pmc_info_list:
            stems.add(pi.File[0:6])
        to_download = [f for f in xml_files if f[0:6] in stems]
        ftp.close()
        return xml_files

    def download_xml_archive(filename, tmp_dir):
        ftp = get_ftp_connection()
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

    def update_text_content(pmc_info_list, xml_file, tmp_dir):
        # Filter PMCIDs to store to only those contained in this XML file
        stem = xml_file[0:6]
        pmc_to_store = [p for p in pmc_info_list if p.File.startswith(stem)]
        # Get the text ref IDs by PMCID
        pmcid_tr_dict = dict(db.get_text_refs_by_pmcid(
                                 tuple([pi.PMCID for pi in pmc_to_store])))
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
                if os.path.exists(xml_path):
                    # Look up the text_ref_id for this PMCID
                    text_ref_id = pmcid_tr_dict[pi.PMCID]
                    # Read the XML file in text mode
                    with open(xml_path, 'rt') as f:
                        content = f.read()
                    # Compress the content
                    content_gz = zip_string(content)
                    # Add to our CSV rows
                    content_block_rows.append([text_ref_id,
                                               b'pmc_auth_xml', content_gz])
                else:
                    print("Could not find file %s" % xml_path)
            t1 = time.time()
            conn = db.get_connection()
            cols = ('text_ref_id', 'content_type', 'content')
            mgr = pgcopy.CopyManager(conn, 'text_content', cols)
            mgr.copy(content_block_rows, BytesIO)

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
    pmc_info_list = get_file_info()
    # 2. Add text_refs to database for any we don't currently have indexed
    update_text_refs(pmc_info_list)
    # 3. Find out which ones are missing auth_xml content in the databse
    pmc_info_missing = get_missing_auth_xml_pmcids(pmc_info_list)
    # 4. Figure out which archives we'll need to download
    xml_list = get_xml_archives_to_download(pmc_info_missing)
    # 5. To save space, download and extract files one at a time
    #xml_list = ['PMC004XXXXXX.xml.tar.gz']
    for xml_file in xml_list:
        tmp_dir = tempfile.mkdtemp(prefix='tmpIndra', dir='.')
        #tmp_dir = 'tmpIndrakgn7bm2j'
        download_xml_archive(xml_file, tmp_dir)
        update_text_content(pmc_info_list, xml_file, tmp_dir)
        shutil.rmtree(tmp_dir)

if __name__ == '__main__':
    #db.drop_tables()
    #db.create_tables()
    initialize_pmc_manuscripts()
    #db.insert_reach('3', '1.3.3', "{'foo': 'bar'}")


"""
tmp_dir = 'tmpIndraipij964n'
filepath = 'tmpIndraipij964n/PMC002XXXXXX.xml.tar.gz'
tf = tarfile.open(name=filepath)
tf.extract(

# Select full list of PMCIDs from database
# Get full list of papers from oa_file_list.csv
# Identify intersection
# Get XMLs for remaining set

# For each file in the list, check the database, and if it's in the DB already
# skip; for the others, obtain.
"""
