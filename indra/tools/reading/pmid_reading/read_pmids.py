from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import argparse
import stat
import random
import sys
import re
import tempfile
import shutil
import subprocess
import glob
import json
import pickle
import functools
import signal
import multiprocessing as mp
from datetime import datetime
from collections import Counter
from platform import system
import logging
from indra import get_config

logger = logging.getLogger('runreader')


parser = argparse.ArgumentParser(
    description=('Apply NLP readers to the content available for a list of '
                 'pmids.')
    )
parser.add_argument(
    '-r', '--reader',
    choices=['reach', 'sparser', 'all'],
    default=['all'],
    dest='readers',
    nargs='+',
    help='Choose which reader(s) to use.'
    )
parser.add_argument(
    '-u', '--upload_json',
    dest='upload_json',
    action='store_true',
    help=('Option to simply upload previously read json files. Overrides -r '
          'option, so no reading will be done.')
    )
parser.add_argument(
    '--force_fulltext',
    dest='force_fulltext',
    action='store_true',
    help='Option to force reading of the full text.'
    )
parser.add_argument(
    '--force_read',
    dest='force_read',
    action='store_true',
    help='Option to force the reader to reread everything.'
    )
parser.add_argument(
    '-n', '--num_cores',
    dest='num_cores',
    default=1,
    type=int,
    help='Select the number of cores you want to use.'
    )
parser.add_argument(
    '-v', '--verbose',
    dest='verbose',
    action='store_true',
    help='Show more output to screen.'
    )
parser.add_argument(
    '-m', '--messy',
    dest='cleanup',
    action='store_false',
    help='Choose to not clean up after run.'
    )
parser.add_argument(
    '-s', '--start_index',
    dest='start_index',
    type=int,
    help='Select the first pmid in the list to start reading.',
    default=0
    )
parser.add_argument(
    '-e', '--end_index',
    dest='end_index',
    type=int,
    help='Select the last pmid in the list to read.',
    default=None
    )
parser.add_argument(
    '--shuffle',
    dest='shuffle',
    action='store_true',
    help=('Select a random sample of the pmids provided. -s/--start_index '
          'will be ingored, and -e/--end_index will set the number of '
          'samples to take.')
    )
parser.add_argument(
    '-o', '--outdir',
    dest='out_dir',
    default=None,
    help=('The output directory where stuff is written. This is only a '
          'temporary directory when reading. By default this will be the'
          '"<basename>_out".')
    )
parser.add_argument(
    dest='basename',
    help='The name of this job.'
    )
parser.add_argument(
    dest='pmid_list_file',
    help=('Path to a file containing a list of line separated pmids for the '
          'articles to be read.')
    )
if __name__ == '__main__':
    args = parser.parse_args()

from indra.sources import reach
from indra.literature import pmc_client, s3_client, get_full_text, \
                             elsevier_client
from indra.sources.sparser import sparser_api as sparser


# Version 1: If JSON is not available, get content and store;
#       assume force_read is False
# Version 1.5: If JSON is not available, get content and store;
#       check for force_read
# Version 2: If JSON is available, return JSON or process
# it and return statements (process it?)


#==============================================================================
# LOADING -- the following are methods to load the content to be read.
#==============================================================================


def join_json_files(prefix):
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


def download_from_s3(pmid, reader='all', input_dir=None, reader_version=None,
                     force_read=False, force_fulltext=False):
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


def get_content_to_read(pmid_list, start_index, end_index, tmp_dir, num_cores,
                        force_fulltext, force_read, reader, reader_version):
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
        download_from_s3,
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


#==============================================================================
# SPARSER -- The following are methods to  process content with sparser.
#==============================================================================


def _timeout_handler(signum, frame):
    raise Exception('Timeout')


def read_pmid(pmid, source, cont_path, sparser_version, outbuf=None,
              cleanup=True):
    "Run sparser on a single pmid."
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(60)
    try:
        if (source is 'content_not_found'
           or source.startswith('unhandled_content_type')
           or source.endswith('failure')):
            logger.info('No content read for %s.' % pmid)
            return  # No real content here.

        if cont_path.endswith('.nxml') and source.startswith('pmc'):
            new_fname = 'PMC%s%d.nxml' % (pmid, mp.current_process().pid)
            os.rename(cont_path, new_fname)

            try:
                sp = sparser.process_nxml_file(
                    new_fname,
                    outbuf=outbuf,
                    cleanup=cleanup
                    )
            finally:
                if cleanup and os.path.exists(new_fname):
                    os.remove(new_fname)
        elif cont_path.endswith('.txt'):
            content_str = ''
            with open(cont_path, 'r') as f:
                content_str = f.read()
            sp = sparser.process_text(
                content_str,
                outbuf=outbuf,
                cleanup=cleanup
                )
        signal.alarm(0)
    except Exception as e:
        logger.error('Failed to process data for %s.' % pmid)
        logger.exception(e)
        signal.alarm(0)
        return

    if sp is None:
        logger.error('Failed to run sparser on pmid: %s.' % pmid)
        return
    # At this point, we rewrite the PMID in the Evidence of Sparser
    # Statements according to the actual PMID that was read.
    sp.set_statements_pmid(pmid)
    s3_client.put_reader_output('sparser', sp.json_stmts, pmid,
                                sparser_version, source)
    return sp.statements


def get_stmts(pmids_unread, cleanup=True, sparser_version=None):
    "Run sparser on the pmids in pmids_unread."
    if sparser_version is None:
        sparser_version = sparser.get_version()
    stmts = {}
    now = datetime.now()
    outbuf_fname = 'sparser_%s_%s.log' % (
        now.strftime('%Y%m%d-%H%M%S'),
        mp.current_process().pid,
        )
    outbuf = open(outbuf_fname, 'wb')
    try:
        for pmid, result in pmids_unread.items():
            logger.info('Reading %s' % pmid)
            source = result['content_source']
            cont_path = result['content_path']
            outbuf.write(('\nReading pmid %s from %s located at %s.\n' % (
                pmid,
                source,
                cont_path
                )).encode('utf-8'))
            outbuf.flush()
            some_stmts = read_pmid(pmid, source, cont_path, sparser_version,
                                   outbuf, cleanup)
            if some_stmts is not None:
                stmts[pmid] = some_stmts
            else:
                continue  # We didn't get any new statements.
    except KeyboardInterrupt as e:
        logger.exception(e)
        logger.info('Caught keyboard interrupt...stopping. \n'
                    'Results so far will be pickled unless '
                    'Keyboard interupt is hit again.')
    finally:
        outbuf.close()
        print("Sparser logs may be found in %s" % outbuf_fname)
    return stmts


def get_stmts_from_cache(pmid):
    json_str = s3_client.get_reader_json_str('sparser', pmid)
    stmts = []
    if json_str is not None:
        stmts = sparser.process_json_dict(json.loads(json_str)).statements
    return {pmid: stmts}


def run_sparser(pmid_list, tmp_dir, num_cores, start_index, end_index,
                force_read, force_fulltext, cleanup=True, verbose=True):
    'Run the sparser reader on the pmids in pmid_list.'
    reader_version = sparser.get_version()
    _, _, _, pmids_read, pmids_unread, _ =\
        get_content_to_read(
            pmid_list, start_index, end_index, tmp_dir, num_cores,
            force_fulltext, force_read, 'sparser', reader_version
            )

    logger.info('Adjusting num cores to length of pmid_list.')
    num_cores = min(len(pmid_list), num_cores)
    logger.info('Adjusted...')
    if num_cores is 1:
        stmts = get_stmts(pmids_unread, cleanup=cleanup)
        stmts.update({pmid: get_stmts_from_cache(pmid)[pmid]
                      for pmid in pmids_read.keys()})
    elif num_cores > 1:
        logger.info("Starting a pool with %d cores." % num_cores)
        pool = mp.Pool(num_cores)
        pmids_to_read = list(pmids_unread.keys())
        N = len(pmids_unread)
        dn = int(N/num_cores)
        logger.info("Breaking pmids into batches.")
        batches = []
        for i in range(num_cores):
            batches.append({
                k: pmids_unread[k]
                for k in pmids_to_read[i*dn:min((i+1)*dn, N)]
                })
        get_stmts_func = functools.partial(
            get_stmts,
            cleanup=cleanup,
            sparser_version=reader_version
            )
        logger.info("Mapping get_stmts onto pool.")
        unread_res = pool.map(get_stmts_func, batches)
        logger.info('len(unread_res)=%d' % len(unread_res))
        read_res = pool.map(get_stmts_from_cache, pmids_read.keys())
        logger.info('len(read_res)=%d' % len(read_res))
        pool.close()
        logger.info('Multiprocessing pool closed.')
        pool.join()
        logger.info('Multiprocessing pool joined.')
        stmts = {
            pmid: stmt_list for res_dict in unread_res + read_res
            for pmid, stmt_list in res_dict.items()
            }
        logger.info('len(stmts)=%d' % len(stmts))

    return (stmts, pmids_unread)


#==============================================================================
# REACH -- The following are methods to process content with reach.
#==============================================================================


REACH_CONF_FMT_FNAME = os.path.join(os.path.dirname(__file__),
                                    '../util/reach_conf_fmt.txt')

REACH_MEM = 5  # GB
MEM_BUFFER = 2  # GB


def process_reach_str(reach_json_str, pmid):
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
    reach_json_str = s3_client.get_reader_json_str('reach', pmid)
    if reach_json_str is None:
        return []
    else:
        return {pmid: process_reach_str(reach_json_str, pmid)}


def upload_reach_readings(pmid, source_type, reader_version, output_dir=None):
    logger.info('Uploading reach result for %s for %s.' % (source_type, pmid))
    # The prefixes should be PMIDs
    prefix_with_path = os.path.join(output_dir, pmid)
    full_json = join_json_files(prefix_with_path)
    # Check that all parts of the JSON could be assembled
    if full_json is None:
        logger.error('REACH output missing JSON for %s' % pmid)
        return {pmid: []}
    # Upload the REACH output to S3
    s3_client.put_reader_output('reach', full_json, pmid, reader_version,
                                source_type)
    return full_json


def upload_process_pmid(pmid_json_tpl):
    pmid, full_json = pmid_json_tpl
    # Process the REACH output with INDRA
    # Convert the JSON object into a string first so that a series of string
    # replacements can happen in the REACH processor
    reach_json_str = json.dumps(full_json)
    return {pmid: process_reach_str(reach_json_str, pmid)}


def upload_process_reach_files(output_dir, pmid_info_dict, reader_version,
                               num_cores):
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
            full_json = upload_reach_readings(
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
    res = pool.map(upload_process_pmid, pmid_json_tuples)
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
    """Run reach on a list of pmids."""
    logger.info('Running REACH with force_read=%s' % force_read)
    logger.info('Running REACH with force_fulltext=%s' % force_fulltext)

    # Get the path to the REACH JAR
    path_to_reach = get_config('REACHPATH')
    if path_to_reach is None or not os.path.exists(path_to_reach):
        logger.warning(
            'Reach path not set or invalid. Check REACHPATH environment var.'
            )
        return {}, {}

    logger.info('Using REACH jar at: %s' % path_to_reach)

    # Get the REACH version
    reach_version = get_config('REACH_VERSION')
    if reach_version is None:
        logger.info('REACH version not set in REACH_VERSION')
        m = re.match('reach-(.*?)\.jar', os.path.basename(path_to_reach))
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
                '-jar', path_to_reach]
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


#==============================================================================
# MAIN -- the main script
#==============================================================================


READER_DICT = {'reach': run_reach, 'sparser': run_sparser}


def get_mem_total():
    if system() == 'Linux':
        with open('/proc/meminfo', 'r') as f:
            lines = f.readlines()
        tot_entry = [line for line in lines if line.startswith('MemTotal')][0]
        ret = int(tot_entry.split(':')[1].replace('kB', '').strip())/10**6
    else:
        ret = None
    return ret


def get_proc_num():
    if system() == 'Linux':
        with open('/proc/cpuinfo', 'r') as f:
            ret = len([
                line for line in f.readlines() if line.startswith('processor')
                ])
    else:
        ret = None
    return ret


def main(args):
    now = datetime.now()    # Set some variables
    if args.upload_json:
        args.readers = 'none'
    out_dir = args.out_dir
    if out_dir is None:
        out_dir = args.basename + '_out'
    made_outdir = False
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        made_outdir = True
    ret = None

    available_cores = get_proc_num()
    if args.num_cores >= available_cores:
        msg = ("You requested %d cores, but only %d available.\n" %
               (args.num_cores, available_cores) +
               "Are you sure you want to proceed? [y/N] > ")
        if sys.version_info.major > 2:
            resp = input(msg)
        else:
            resp = raw_input(msg)
        if resp.lower() not in ['y', 'yes']:
            logger.info("Aborting...")
            return
    try:
        # Option -u: just upload previously read JSON files
        if args.upload_json:
            with open(args.pmid_list_file, 'rb') as f:
                text_sources = pickle.load(f)
            stmts = upload_process_reach_files(
                args.out_dir,
                text_sources,
                args.num_cores
                )
            pickle_file = '%s_stmts.pkl' % args.basename
            with open(pickle_file, 'wb') as f:
                pickle.dump(stmts, f, protocol=2)
            sys.exit()

        # Option -r <reader>: actually read the content.

        # Load the list of PMIDs from the given file
        with open(args.pmid_list_file) as f:
            pmid_list = [line.strip('\n') for line in f.readlines()]

        if args.end_index is None:
            args.end_index = len(pmid_list)

        if args.shuffle:
            pmid_list = random.sample(pmid_list, args.end_index)

        # Do the reading
        readers = []
        if 'all' in args.readers:
            readers = list(READER_DICT.keys())
        else:
            readers = args.readers[:]

        stmts = {}
        for reader in readers:
            run_reader = READER_DICT[reader]
            some_stmts, _ = run_reader(
                pmid_list,
                out_dir,
                args.num_cores,
                args.start_index,
                args.end_index,
                args.force_read,
                args.force_fulltext,
                cleanup=args.cleanup,
                verbose=args.verbose
                )
            stmts[reader] = some_stmts

        N_tot = sum([
            len(stmts[reader][pmid]) for reader in readers
            for pmid in stmts[reader].keys()
            ])
        logger.info('Collected a total of %s statements.' % N_tot)
        for reader in readers:
            N = sum([len(l) for l in stmts[reader].values()])
            logger.info(
                '%s accumulated %d statements.' % (reader.capitalize(), N)
                )

        # Pickle the statements
        if args.end_index is None:
            args.end_index = 'end'
        pickle_file = '%s_stmts_%s-%s.pkl' % (
            args.basename,
            args.start_index,
            args.end_index
            )
        with open(pickle_file, 'wb') as f:
            pickle.dump(stmts, f, protocol=2)
        ret = pickle_file
    finally:
        time_taken = datetime.now() - now
        print('This run took', time_taken)
        time_file = os.path.join(os.path.dirname(__file__), 'time_data.txt')
        with open(time_file, 'a') as f:
            f.write('Started run at %s with args %s lasting %s.\n' %
                    (now, str(args), time_taken))
        if made_outdir and args.cleanup:
            shutil.rmtree(out_dir)
    return ret


if __name__ == '__main__':
    main(args)
