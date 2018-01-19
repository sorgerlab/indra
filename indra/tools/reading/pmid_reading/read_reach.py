"""Methods to process content with REACH."""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import os
import json
import glob
import shutil
import logging
import subprocess
import multiprocessing as mp
from indra.sources import reach
from indra.literature import s3_client
from .util import get_mem_total
from .get_content import get_content_to_read

logger = logging.getLogger('pmid_reading/read_reach')


REACH_CONF_FMT_FNAME = os.path.join(os.path.dirname(__file__),
                                    'reach_conf_fmt.txt')

REACH_MEM = 5  # GB
MEM_BUFFER = 2  # GB


def process_reach_str(reach_json_str, pmid):
    """Process a REACH output JSON string and return Statements."""
    if reach_json_str is None:
        raise ValueError('reach_json_str cannot be None')
    # Run the REACH processor on the JSON
    try:
        reach_proc = reach.process_json_str(reach_json_str, citation=pmid)
    # If there's a problem, skip it
    except Exception as e:
        print("Exception processing %s" % pmid)
        print(e)
        return []
    return reach_proc.statements


def process_reach_from_s3(pmid):
    """Return the given PMID and Statements processed from it.

    This function gets the REACH output from S3, processes it
    and returns a dict with the given PMID and the Statements
    obtained by processing the reading result.
    """
    reach_json_str = s3_client.get_reader_json_str('reach', pmid)
    if reach_json_str is None:
        return []
    else:
        return {pmid: process_reach_str(reach_json_str, pmid)}


def upload_reach_json(pmid, source_type, reader_version, output_dir=None):
    """Join and upload the output of REACH for a given PMID."""
    logger.info('Uploading reach result for %s for %s.' % (source_type, pmid))
    # The prefixes should be PMIDs
    prefix_with_path = os.path.join(output_dir, pmid)
    full_json = join_reach_json_files(prefix_with_path)
    # Check that all parts of the JSON could be assembled
    if full_json is None:
        logger.error('REACH output missing JSON for %s' % pmid)
        return {pmid: []}
    # Upload the REACH output to S3
    s3_client.put_reader_output('reach', full_json, pmid, reader_version,
                                source_type)
    return full_json


def process_reach_output_pmid(pmid_json_tpl):
    """Return a dict of a given PMID with extracted Statements.

    This function processes the JSON output of REACH for a single
    PMID.
    """
    pmid, full_json = pmid_json_tpl
    # Process the REACH output with INDRA
    # Convert the JSON object into a string first so that a series of string
    # replacements can happen in the REACH processor
    reach_json_str = json.dumps(full_json)
    return {pmid: process_reach_str(reach_json_str, pmid)}


def upload_process_reach_files(output_dir, pmid_info_dict, reader_version,
                               num_cores):
    """Upload and process all reading output from REACH."""
    # At this point, we have a directory full of JSON files
    # Collect all the prefixes into a set, then iterate over the prefixes

    # Collect prefixes
    json_files = glob.glob(os.path.join(output_dir, '*.json'))
    json_prefixes = set([])
    for json_file in json_files:
        filename = os.path.basename(json_file)
        prefix = filename.split('.')[0]
        json_prefixes.add(prefix)
    # Make a list with PMID and source_text info
    logger.info("Uploading reading results for reach.")
    pmid_json_tuples = []
    for json_prefix in json_prefixes:
        try:
            full_json = upload_reach_json(
                json_prefix,
                pmid_info_dict[json_prefix].get('content_source'),
                reader_version,
                output_dir
                )
            pmid_json_tuples.append((json_prefix, full_json))
        except Exception as e:
            logger.error("Caught an exception while trying to upload reach "
                         "reading results onto s3 for %s." % json_prefix)
            logger.exception(e)
    # Create a multiprocessing pool
    logger.info('Creating a multiprocessing pool with %d cores' % num_cores)
    # Get a multiprocessing pool.
    pool = mp.Pool(num_cores)
    logger.info('Processing local REACH JSON files')
    res = pool.map(process_reach_output_pmid, pmid_json_tuples)
    stmts_by_pmid = {
        pmid: stmts for res_dict in res for pmid, stmts in res_dict.items()
        }
    pool.close()
    logger.info('Multiprocessing pool closed.')
    pool.join()
    logger.info('Multiprocessing pool joined.')
    """
    logger.info('Uploaded REACH JSON for %d files to S3 (%d failures)' %
        (num_uploaded, num_failures))
    failures_file = os.path.join(output_dir, 'failures.txt')
    with open(failures_file, 'wt') as f:
        for fail in failures:
            f.write('%s\n' % fail)
    """
    return stmts_by_pmid


def run_reach(pmid_list, base_dir, num_cores, start_index, end_index,
              force_read, force_fulltext, cleanup=False, verbose=True):
    """Run reach on a list of PMIDs.

    As opposed to Sparser, REACH is called only once to run on all the
    PMIDs.
    """
    logger.info('Running REACH with force_read=%s' % force_read)
    logger.info('Running REACH with force_fulltext=%s' % force_fulltext)

    # Get the path to the reach directory.
    path_to_reach = os.environ.get('REACHPATH', None)
    if path_to_reach is None or not os.path.exists(path_to_reach):
        logger.warning(
            'Reach path not set or invalid. Check REACHPATH environment var.'
            )
        return {}, {}
    patt = re.compile('reach-(.*?)\.jar')

    # Find the jar file.
    for fname in os.listdir(path_to_reach):
        m = patt.match(fname)
        if m is not None:
            reach_ex = os.path.join(path_to_reach, fname)
            break
    else:
        logger.warning("Could not find reach jar in reach dir.")
        return {}, {}

    logger.info('Using REACH jar at: %s' % reach_ex)

    # Get the reach version.
    reach_version = os.environ.get('REACH_VERSION', None)
    if reach_version is None:
        logger.info('REACH version not set in REACH_VERSION')
        reach_version = re.sub('-SNAP.*?$', '', m.groups()[0])

    logger.info('Using REACH version: %s' % reach_version)

    tmp_dir, _, output_dir, pmids_read, pmids_unread, num_found =\
        get_content_to_read(
            pmid_list, start_index, end_index, base_dir, num_cores,
            force_fulltext, force_read, 'reach', reach_version
            )

    stmts = {}
    mem_tot = get_mem_total()
    if mem_tot is not None and mem_tot <= REACH_MEM + MEM_BUFFER:
        logger.error(
            "Too little memory to run reach. At least %s required." %
            REACH_MEM + MEM_BUFFER
            )
        logger.info("REACH not run.")
    elif len(pmids_unread) > 0 and num_found > 0:
        # Create the REACH configuration file
        with open(REACH_CONF_FMT_FNAME, 'r') as fmt_file:
            conf_file_path = os.path.join(tmp_dir, 'indra.conf')
            with open(conf_file_path, 'w') as conf_file:
                conf_file.write(
                    fmt_file.read().format(tmp_dir=os.path.abspath(tmp_dir),
                                           num_cores=num_cores,
                                           loglevel='INFO')
                    )

        # Run REACH!
        logger.info("Beginning reach.")
        args = ['java', '-Xmx24000m', '-Dconfig.file=%s' % conf_file_path,
                '-jar', reach_ex]
        p = subprocess.Popen(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if verbose:
            for line in iter(p.stdout.readline, b''):
                logger.info(line)
        p_out, p_err = p.communicate()
        if p.returncode:
            logger.error('Problem running REACH:')
            logger.error('Stdout: %s' % p_out.decode('utf-8'))
            logger.error('Stderr: %s' % p_err.decode('utf-8'))
            raise Exception('REACH crashed')

        # Process JSON files from local file system, process to INDRA
        # Statements and upload to S3
        some_stmts = upload_process_reach_files(
            output_dir,
            pmids_unread,
            reach_version,
            num_cores
            )
        stmts.update(some_stmts)
        # Delete the tmp directory if desired
        if cleanup:
            shutil.rmtree(tmp_dir)

    # Create a new multiprocessing pool for processing the REACH JSON
    # files previously cached on S3
    logger.info('Creating multiprocessing pool with %d cpus' % num_cores)
    pool = mp.Pool(num_cores)

    # Download and process the JSON files on S3
    logger.info('Processing REACH JSON from S3 in parallel')
    res = pool.map(process_reach_from_s3, pmids_read.keys())
    pool.close()
    logger.info('Multiprocessing pool closed.')
    pool.join()
    logger.info('Multiprocessing pool joined.')
    s3_stmts = {
        pmid: stmt_list for res_dict in res
        for pmid, stmt_list in res_dict.items()
        }
    stmts.update(s3_stmts)

    # Save the list of PMIDs with no content found on S3/literature client
    '''
    content_not_found_file = os.path.join(tmp_dir, 'content_not_found.txt')
    with open(content_not_found_file, 'wt') as f:
        for c in content_not_found:
            f.write('%s\n' % c)
    '''
    return stmts, pmids_unread


def join_reach_json_files(prefix):
    """Join different REACH output JSON files into a single JSON object.

    The output of REACH is broken into three files that need to be joined
    before processing. Specifically, there will be three files of the form:
    `<prefix>.uaz.<subcategory>.json`.

    Parameters
    ----------
    prefix : str
        The absolute path up to the extensions that reach will add.

    Returns
    -------
    json_obj : dict
        The result of joining the files, keyed by the three subcategories.
    """
    try:
        with open(prefix + '.uaz.entities.json', 'rt') as f:
            entities = json.load(f)
        with open(prefix + '.uaz.events.json', 'rt') as f:
            events = json.load(f)
        with open(prefix + '.uaz.sentences.json', 'rt') as f:
            sentences = json.load(f)
    except IOError as e:
        logger.error(
            'Failed to open JSON files for %s; REACH error?' % prefix
            )
        logger.exception(e)
        return None
    return {'events': events, 'entities': entities, 'sentences': sentences}

