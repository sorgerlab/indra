import os
import shutil
import tarfile
import tempfile
from ftplib import FTP
from collections import namedtuple
import db_api

pmc_ftp_url = 'ftp.ncbi.nlm.nih.gov'
blocksize=33554432 # Chunk size recommended by NCBI

def initialize_pmc_manuscripts():

    auth_dir = '/pub/pmc/manuscript'
    #temp_dir = tempfile.mkdtemp(prefix='tmpIndra', dir='.')
    temp_dir = 'tmpIndrari7qqk_i'
    # Get an FTP connection
    ftp = FTP(pmc_ftp_url)
    ftp.login()
    # Change to the manuscripts directory
    ftp.cwd(auth_dir)
    # Some variables for meaningful progress messages
    stored_bytes = 0
    log_pcts = list(range(0, 101, 5))
    # Get the list of files from the CSV file
    filelist_bytes = []
    print("Downloading filelist.csv")
    ftp.retrbinary('RETR filelist.csv',
                  callback=lambda b: filelist_bytes.append(b),
                  blocksize=blocksize)
    filelist_csv = b''.join(filelist_bytes).decode('ascii').split('\n')
    # Namedtuple for working with PMC info entries
    print("Processing filelist.csv")
    PmcInfo = namedtuple('PmcInfo', ('File', 'PMCID', 'PMID', 'MID'))
    # Process the file info (skip the header line)
    pmc_info_list = [PmcInfo(*line.split(',')) for line in filelist_csv
                     if line][1:]
    """
    # Function to write to local file with progress updates
    def write_to_file(fp, b, total_size):
        nonlocal stored_bytes, log_pcts
        fp.write(b)
        stored_bytes += len(b)
        pct_complete = round(100 * (stored_bytes / float(total_size)))
        if pct_complete in log_pcts:
            print('%s: %s%% complete' % (filename, pct_complete))
            log_pcts.remove(pct_complete)

    filename = 'PMC002XXXXXX.xml.tar.gz'
    outfilepath = os.path.join(temp_dir, filename)
    filesize = ftp.size(filename)
    print("Getting %s" % filename)
    with open(outfilepath, 'wb') as f:
        ftp.retrbinary('RETR %s' % filename,
                       callback=lambda b: write_to_file(f, b, filesize),
                       blocksize=blocksize)
    ftp.close()
    # Extract all files in the TAR archive
    tf = tarfile.open(outfilepath)
    print("Extracting all files from %s" % outfilepath)
    tf.extractall(path=tmp_dir)
    # Now that we've extracted everything, walk the directory tree and upload
    # content to the INDRA database
    """
    counter = 0
    for pmc_info in pmc_info_list:
        if counter > 5:
            break
        xml_path = os.path.join(temp_dir, pmc_info.File)
        if os.path.exists(xml_path):
            print("Found file %s" % xml_path)
            db_api.add_text_ref(source='pmc', pmid=pmc_info.PMID,
                                pmcid=pmc_info.PMCID, manuscript_id=pmc_info.MID)
            # Open in text mode
            with open(xml_path, 'rt') as f:
                content = f.read()
            db_api.add_text_content_by_pmid(pmc_info.PMID, 'pmc_auth_xml',
                                            content)
        counter += 1

def update_pmc_manuscripts():
    pass



"""


temp_dir = 'tmpIndraipij964n'
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
