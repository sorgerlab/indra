"""Methods to load literature content to be read."""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import stat
import pickle
import logging
import tempfile
import functools
import multiprocessing as mp
from collections import Counter
from indra.literature import s3_client, elsevier_client


logger = logging.getLogger('pmid_reading/get_content')


def get_content_to_read(pmid_list, start_index, end_index, tmp_dir, num_cores,
                        force_fulltext, force_read, reader, reader_version):
    """Download and return information on text content to be read.

    This function takes a list of PMIDs and uses the download_from_s3 function
    to download cached literature content from S3. If any new content download
    is needed, it is handled within download_from_s3.
    """
    if end_index is None or end_index > len(pmid_list):
        end_index = len(pmid_list)
    pmids_in_range = pmid_list[start_index:end_index]

    # Create the temp directories for input and output
    base_dir = tempfile.mkdtemp(prefix='read_%s_to_%s_' %
                                (start_index, end_index), dir=tmp_dir)
    # Make the temp directory writeable by REACH
    os.chmod(base_dir, stat.S_IRWXO | stat.S_IRWXU | stat.S_IRWXG)
    input_dir = os.path.join(base_dir, 'input')
    output_dir = os.path.join(base_dir, 'output')
    os.makedirs(input_dir)
    os.makedirs(output_dir)

    download_from_s3_func = functools.partial(
        download_content_from_s3,
        input_dir=input_dir,
        reader=reader,
        reader_version=reader_version,
        force_read=force_read,
        force_fulltext=force_fulltext
        )

    if num_cores > 1:
        # Get content using a multiprocessing pool
        logger.info('Creating multiprocessing pool with %d cpus' % num_cores)
        pool = mp.Pool(num_cores)
        logger.info('Getting content for PMIDs in parallel')
        res = pool.map(download_from_s3_func, pmids_in_range)
        pool.close()  # Wait for procs to end.
        logger.info('Multiprocessing pool closed.')
        pool.join()
        logger.info('Multiprocessing pool joined.')
    else:
        res = list(map(download_from_s3_func, pmids_in_range))

    # Combine the results into a single dict
    pmid_results = {
        pmid: results for pmid_dict in res
        for pmid, results in pmid_dict.items()
        }
    # Tabulate and log content results here
    pmids_read = {
        pmid: result for pmid, result in pmid_results.items()
        if result.get('reader_version') == reader_version
        }
    pmids_unread = {
        pmid: pmid_results[pmid]
        for pmid in set(pmid_results.keys()).difference(set(pmids_read.keys()))
        }
    logger.info(
        '%d / %d papers already read with %s %s' %
        (len(pmids_read), len(pmid_results), reader, reader_version)
        )
    num_found = len([
        pmid for pmid in pmids_unread
        if pmids_unread[pmid].get('content_path') is not None
        ])
    logger.info(
        'Retrieved content for %d / %d papers to be read' %
        (num_found, len(pmids_unread))
        )

    # Tabulate sources and log in sorted order
    content_source_list = [
        pmid_dict.get('content_source') for pmid_dict in pmids_unread.values()
        ]
    content_source_counter = Counter(content_source_list)
    content_source_list = [
        (source, count) for source, count in content_source_counter.items()
        ]
    content_source_list.sort(key=lambda x: x[1], reverse=True)
    if content_source_list:
        logger.info('Content sources:')
        for source, count in content_source_list:
            logger.info('%s: %d' % (source, count))

    # Save text sources
    logger.info('Saving text sources...')
    text_source_file = os.path.join(base_dir, 'content_types.pkl')
    with open(text_source_file, 'wb') as f:
        pickle.dump(pmids_unread, f, protocol=2)
    logger.info('Text sources saved.')

    return base_dir, input_dir, output_dir, pmids_read, pmids_unread, num_found


def download_content_from_s3(pmid, reader='all', input_dir=None,
                             reader_version=None, force_read=False,
                             force_fulltext=False):
    """Return content and metadata for content for a single PMID.

    This function attempts to get cached literature content from S3, and if
    not available."""
    logger.info(('Downloading %s from S3, force_read=%s, force_fulltext=%s '
                 'reader_version=%s') % (pmid, force_read, force_fulltext,
                                         reader_version))
    if input_dir is None:
        raise ValueError('input_dir must be defined')

    # First define the text retrieval function
    def get_text():
        # full_pmid = s3_client.check_pmid(pmid)
        # Look for the full text
        content, content_type = s3_client.get_upload_content(
            pmid,
            force_fulltext_lookup=force_fulltext
            )
        content_path = None
        # Write the contents to a file
        if content_type is None or content is None:
            # No content found on S3, skipping
            content_source = 'content_not_found'
        elif content_type == 'pmc_oa_xml':
            content_source = 'pmc_oa_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_auth_xml':
            content_source = 'pmc_auth_xml'
            content_path = os.path.join(input_dir, '%s.nxml' % pmid)
        elif content_type == 'pmc_oa_txt':
            content_source = 'pmc_oa_txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'elsevier_xml':
            content = elsevier_client.extract_text(content)
            # Couldn't get text from Elsevier XML
            if content is None:
                content_source = 'elsevier_extract_text_failure'
            else:
                content_source = 'elsevier_xml'
                content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'txt':
            content_source = 'txt'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        elif content_type == 'abstract':
            content_source = 'abstract'
            content_path = os.path.join(input_dir, '%s.txt' % pmid)
        # Unhandled content type, skipping
        else:
            content_source = 'unhandled_content_type_%s' % content_type
        # If we got content, write the content to a file with the appropriate
        # extension
        if content_path:
            with open(content_path, 'wb') as f:
                # The XML string is Unicode
                enc = content.encode('utf-8')
                f.write(enc)
        # Return dict of results for this PMID
        result = {pmid: {'content_source': content_source,
                         'content_path': content_path}}
        return result

    # If we're forcing a read regardless of whether there is cached REACH
    # output, then we download the text content
    if force_read or reader_version is None:
        return get_text()
    # If not, look for REACH JSON on S3
    reader_version_s3, read_source_text = \
        s3_client.get_reader_metadata(reader, pmid)
    # Found it, same version, no need to get text
    if (reader_version_s3 is not None
       and reader_version_s3 == reader_version):
        result = {pmid: {
            'reader_version': reader_version_s3,
            'reach_source_text': read_source_text
            }}
    # Found it, different version, get the text
    else:
        result = get_text()
        result[pmid].update({'reader_version': reader_version_s3,
                             'reach_source_text': read_source_text})
    return result
