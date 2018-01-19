"""Methods to  process content with Sparser."""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import signal
import logging
import functools
import multiprocessing as mp
from datetime import datetime
from indra.sources import sparser
from indra.literature import s3_client
from .get_content import get_content_to_read

logger = logging.getLogger('pmid_reading/read_sparser')


def run_sparser(pmid_list, tmp_dir, num_cores, start_index, end_index,
                force_read, force_fulltext, cleanup=True, verbose=True):
    """Run the Sparser reader on a given list of PMIDs."""
    # Get the version of the Sparser reader
    reader_version = sparser.get_version()

    # Check which PMIDs have already been read and which haven't yet
    _, _, _, pmids_read, pmids_unread, _ = \
        get_content_to_read(
            pmid_list, start_index, end_index, tmp_dir, num_cores,
            force_fulltext, force_read, 'sparser', reader_version
            )

    # Adjust number of cores to do reading with
    logger.info('Adjusting num cores to length of pmids_unread.')
    num_cores = min(len(pmids_unread), num_cores)
    logger.info('Adjusted to use %d cores...' % num_cores)

    # If it's a single core, just call read_pmids_sparser directly
    if num_cores == 1:
        # Get Statements by reading
        stmts_from_reading = read_pmids(pmids_unread, cleanup=cleanup)
        # Get Statements from cache
        stmts_from_cache = {pmid: process_from_s3(pmid)[pmid]
                            for pmid in pmids_read.keys()}
        # Combine the two dicts
        stmts = stmts_from_reading
        stmts.update(stmts_from_cache)
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
        read_pmids_func = functools.partial(
            read_pmids,
            cleanup=cleanup,
            sparser_version=reader_version
            )
        logger.info("Mapping read_pmids_sparser onto pool.")
        # Get results of reading as dict of PMIDs with Statements
        unread_res = pool.map(read_pmids_func, batches)
        # Get results from cache as dict of PMIDs with Statements
        read_res = pool.map(process_from_s3, pmids_read.keys())
        pool.close()
        logger.info('Multiprocessing pool closed.')
        pool.join()
        logger.info('Multiprocessing pool joined.')
        # Combine the reading and cache results into a single dict
        # of PMIDs with Statements
        stmts = {
            pmid: stmt_list for res_dict in unread_res + read_res
            for pmid, stmt_list in res_dict.items()
            }

    return (stmts, pmids_unread)


def read_pmids(pmids, cleanup=True, sparser_version=None):
    """Run sparser on a list of PMIDs and return a dict of Statements

    Given a list of PMIDs, this function calls read_one_pmid
    to read each PMID, and builds up a dict of Statements with
    the PMIDs as keys.
    """
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
        for pmid, result in pmids.items():
            logger.info('Reading %s' % pmid)
            source = result['content_source']
            cont_path = result['content_path']
            outbuf.write(('\nReading pmid %s from %s located at %s.\n' % (
                pmid,
                source,
                cont_path
                )).encode('utf-8'))
            outbuf.flush()
            some_stmts = read_one_pmid(pmid, source, cont_path,
                                               sparser_version,
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


def _timeout_handler(signum, frame):
    raise Exception('Timeout')


def read_one_pmid(pmid, source, cont_path, sparser_version, outbuf=None,
                  cleanup=True):
    """Run Sparser on a single PMID and return a list of Statements.

    This function runs Sparser on a single paper and caches the result of
    reading (a JSON file) on S3.
    """
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
    s3_client.put_reader_output('sparser', sp.json_stmts, pmid,
                                sparser_version, source)
    return sp.statements


def process_from_s3(pmid):
    """Return Statements that Sparser extracted for the given PMID.

    The result is returned as a dict with the given PMID as the key
    and a list of Statements as the value.
    """
    json_str = s3_client.get_reader_json_str('sparser', pmid)
    stmts = []
    if json_str is not None:
        sp = sparser.process_json_dict(json.loads(json_str))
        if sp is not None:
            stmts = sp.statements
    return {pmid: stmts}
