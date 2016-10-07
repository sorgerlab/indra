from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import botocore
from indra.literature import s3_client
from indra.literature import get_full_text
from indra.literature import pubmed_client
import logging

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: %s pmid_list" % sys.argv[0])
        sys.exit()

    pmid_list = sys.argv[1]
    logger = logging.getLogger('upload_content')

    with open(pmid_list) as f:
        pmids = [line.strip('\n') for line in f.readlines()]

    # Search for full text and abstract on S3 and store info about what needs
    # to be uploaded (fulltext, abstract) in a list of tuples.
    pmids_to_upload = []
    logger.info("-- Checking for %d full texts --" % len(pmids))
    for ix, pmid in enumerate(pmids[0:10]):
        s3_client.get_upload_content(pmid)
