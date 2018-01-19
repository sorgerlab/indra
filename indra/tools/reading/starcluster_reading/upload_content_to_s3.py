from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import sys
import botocore
from indra.literature import s3_client
from indra.literature import get_full_text
from indra.literature import pubmed_client
import logging

if __name__ == '__main__':

    usage = "Usage: %s pmid_list start_index end_index [force_fulltext]" % \
             sys.argv[0]

    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print(usage)
        sys.exit()
    if len(sys.argv) == 5 and sys.argv[4] != 'force_fulltext':
        print(usage)
        sys.exit()
    elif len(sys.argv) == 5:
        force_fulltext = True
    else:
        force_fulltext = False

    logger = logging.getLogger('upload_content')

    pmid_list = sys.argv[1]
    start_index = int(sys.argv[2])
    end_index = int(sys.argv[3])

    with open(pmid_list) as f:
        pmids = [line.strip('\n') for line in f.readlines()]

    # Search for full text and abstract on S3 and store info about what needs
    # to be uploaded (fulltext, abstract) in a list of tuples.
    pmids_to_upload = []
    logger.info("-- Checking for %d full texts --" % len(pmids))
    if end_index > len(pmids):
        end_index = len(pmids)
    for ix, pmid in enumerate(pmids[start_index:end_index]):
        logger.info("--- %d: %s ---" % ((start_index + ix), pmid))
        s3_client.get_upload_content(pmid,
                                     force_fulltext_lookup=force_fulltext)
