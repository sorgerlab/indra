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
    for ix, pmid in enumerate(pmids):
        logger.info("%d: %s" % (ix + 1, pmid))
        upload_info = [pmid]
        # Check for any full text. If there is none, try to get and upload.
        ft_prefix = s3_client.get_pmid_key(pmid) + '/fulltext/'
        ft_objs = s3_client.filter_keys(ft_prefix)
        if len(ft_objs) == 0:
            upload_info.append('fulltext')

        # Check if we have an abstract
        abs_key = s3_client.get_pmid_key(pmid) + '/abstract'
        if not s3_client.check_key(abs_key):
            upload_info.append('abstract')
        # Check if we need to upload anything
        if len(upload_info) > 1:
            pmids_to_upload.append(upload_info)
        with open('upload_info.csv', 'a') as f:
            f.write(','.join(upload_info) + '\n')

    # Now, iterate over pmids and upload content to S3
    logger.info("\n-- Uploading content for %d papers --" %
                len(pmids_to_upload))
    for ix, upload_info in enumerate(pmids_to_upload):
        logger.info("%d: %s" % (ix + 1, pmid))
        pmid = upload_info[0]
        abstract = None
        if 'fulltext' in upload_info:
            (content, content_type) = get_full_text(pmid, idtype='pmid')
            # If content is None that means that the ID lookup failed in some
            # way so we will skip this paper
            if content is None:
                continue

            if content_type == 'nxml':
                logger.info("Uploading nxml for %s" % pmid)
                s3_client.put_full_text(pmid, content,
                                        full_text_type='pmc_oa_xml')
            elif content_type == 'txt':
                logger.info("Uploading txt for %s" % pmid)
                s3_client.put_full_text(pmid, content,
                                        full_text_type='txt')
            elif content_type == 'abstract':
                abstract = content
        # Now check if we need to upload the abstract (which we may have already
        # gotten)
        if 'abstract' in upload_info:
            # Check to see if we already got the abstract
            if abstract is not None:
                logger.info("Putting abstract for %s" % pmid)
                s3_client.put_abstract(pmid, abstract)
            else:
                abstract = pubmed_client.get_abstract(pmid)
                if abstract is not None:
                    logger.info("Putting abstract for %s" % pmid)
                    s3_client.put_abstract(pmid, abstract)
        with open('pmids_done.txt', 'a') as f:
            f.write('%s\n' % pmid)

