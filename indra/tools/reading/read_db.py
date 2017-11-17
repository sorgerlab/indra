"""This module provides essential tools to run reading using indra's own
database. This may also be run as a script; for details run:
`python read_pmids_db --help`
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import logging
import random
import zlib
import pickle
from docutils.io import InputError
from math import log10, floor
from os.path import join as pjoin

from indra.db import get_primary_db, formats, texttypes
from indra.db import sql_expressions as sql

from indra.tools.reading.readers import get_readers, ReadingData, _get_dir
from indra.tools.reading.script_tools import get_parser, make_statements, \
                                             StatementData

logger = logging.getLogger('make_db_readings')
if __name__ == '__main__':
    parser = get_parser(
        'A tool to read and process content from the database.',
        'A file containing a list of ids of the form <id_type>:<id>.'
        )
    parser.add_argument(
        '-p', '--pickle',
        help='Pickle all results and save in .pkl files.',
        action='store_true'
        )
    parser.add_argument(
        '--no_reading_upload',
        help='Choose not to upload the reading output to the database.',
        action='store_true'
        )
    parser.add_argument(
        '--no_statement_upload',
        help='Choose not to upload the statements to the databse.',
        action='store_true'
        )
    parser.add_argument(
        '--force_read',
        help=('Make the reader read all the content. Otherwise, the system '
              'will simply process existing readings where possible.'),
        action='store_true'
        )
    parser.add_argument(
        '--force_fulltext',
        help='Make the reader only read full text from the database.',
        action='store_true'
        )
    parser.add_argument(
        '--no_read',
        help='Only create statements using existing reading results.',
        action='store_true'
        )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)


# =============================================================================
# Useful functions
# =============================================================================


def _convert_id_entry(id_entry, allowed_types=None):
    ret = [s.strip() for s in id_entry.split(':')]
    if len(ret) != 2:
        raise InputError('Improperly formatted id entry: \"%s\"' % id_entry)
    ret[0] = ret[0].lower()
    if allowed_types is not None and ret[0] not in allowed_types:
        raise InputError('Invalid id type: \"%s\"' % ret[0])
    return ret


def _enrich_reading_data(reading_data_iter, db=None):
    """Get db ids for all ReadingData objects that correspond to a db ref.

    Note that the objects are modified IN PLACE, so nothing is returned, and if
    a copy of the objects is passed as an argument, this function will have no
    effect. This does nothing if the readings are not in the database.
    """
    logger.debug("Enriching the reading data with database refs.")
    if db is None:
        db = get_primary_db()
    possible_matches = db.select_all(
        'readings',
        db.Readings.text_content_id.in_([rd.tcid for rd in reading_data_iter
                                         if rd.reading_id is None])
        )
    for rdata in reading_data_iter:
        for reading in possible_matches:
            if rdata.matches(reading):
                rdata.reading_id = reading.id
                break
    return


# =============================================================================
# Content Retrieval
# =============================================================================


def get_clauses(id_str_list, db):
    """Get a list of clauses to be passed to a db query."""
    id_types = db.TextRef.__table__.columns.keys()
    id_dict = {id_type: [] for id_type in id_types}
    for id_entry in id_str_list:
        id_type, id_val = _convert_id_entry(id_entry, id_types)
        id_dict[id_type].append(id_val)
    id_condition_list = [getattr(db.TextRef, id_type).in_(id_list)
                         for id_type, id_list in id_dict.items()
                         if len(id_list)]
    return [sql.or_(*id_condition_list)]


def get_text_content_summary_string(q, db, num_ids=None):
    """Create a table with some summary data for a query."""
    N_tot = q.count()
    if num_ids is not None:
        logger.info("Found %d text content entires out of %d ids."
                    % (N_tot, num_ids))
    if N_tot > 0:
        log_n = floor(log10(N_tot))
    else:
        log_n = 1
    cols = list(formats.values()) + ['tot']
    col_fmt = ' %%%ds' % max(4, log_n)
    cols_strs = [col_fmt % fmt for fmt in cols]
    ret_str = 'Summary Statisitics:\n' + ' '*10 + ''.join(cols_strs) + '\n'
    col_counts = dict.fromkeys(formats.values())
    col_counts['tot'] = []
    for texttype in texttypes.values():
        line = '%8s: ' % texttype
        counts = []
        for text_format in cols[:-1]:
            if col_counts[text_format] is None:
                col_counts[text_format] = []
            c = q.filter(
                db.TextContent.text_type == texttype,
                db.TextContent.format == text_format
                ).count()
            line += col_fmt % c
            counts.append(c)
            col_counts[text_format].append(c)
        line += col_fmt % sum(counts)
        ret_str += line + '\n'
        col_counts['tot'].append(sum(counts))
    ret_str += '%8s: ' % 'total' + ''.join([col_fmt % sum(col_counts[col])
                                            for col in cols])
    return ret_str


def get_content_query(id_str_list, readers, db=None, force_fulltext=False,
                      force_read=False):
    """Load all the content that will be read."""
    if db is None:
        db = get_primary_db()
    logger.debug("Got db handle.")

    tc_tr_binding = db.TextContent.text_ref_id == db.TextRef.id
    rd_tc_binding = db.Readings.text_content_id == db.TextContent.id
    general_clauses = [tc_tr_binding]
    if force_fulltext:
        general_clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)
    id_clauses = get_clauses(id_str_list, db)
    logger.debug("Generated %d id clauses." % len(id_clauses))

    if id_clauses:
        # Get the text content query object
        tc_query = db.filter_query(
            db.TextContent,
            *(general_clauses + id_clauses)
            ).distinct()
        try:
            logger.debug("Going to try to make a nice summary...")
            logger.info(get_text_content_summary_string(tc_query, db))
        except Exception:
            logger.info("Could not print summary of results.")

        if not force_read:
            logger.debug("Getting content to be read.")
            # Each sub query is a set of content that has been read by one of
            # the readers.
            tc_q_subs = [tc_query.filter(rd_tc_binding, r.matches_clause(db))
                         for r in readers]
            tc_tbr_query = tc_query.except_(sql.intersect(*tc_q_subs))
        else:
            logger.debug('All content will be read (force_read).')
            tc_tbr_query = tc_query
    else:
        logger.info("Did not retreive content or readings from database.")
        tc_tbr_query = None

    return tc_tbr_query.distinct()


def get_readings_query(id_str_list, readers, db=None, force_fulltext=False):
    """Get all the readings available for the id's in id_str_list."""
    if db is None:
        db = get_primary_db()
    general_clauses = [
        # Bind conditions on readings to conditions on content.
        db.Readings.text_content_id == db.TextContent.id,

        # Bind text content to text refs
        db.TextContent.text_ref_id == db.TextRef.id,

        # Check if at least one of the readers has read the content
        sql.or_(*[reader.matches_clause(db) for reader in readers])
        ]
    if force_fulltext:
        general_clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)
    id_clauses = get_clauses(id_str_list, db)
    if id_clauses:
        readings_query = db.filter_query(
            db.Readings,

            # Bind conditions on readings to conditions on content.
            db.Readings.text_content_id == db.TextContent.id,

            # Bind text content to text refs
            db.TextContent.text_ref_id == db.TextRef.id,

            # Check if at least one of the readers has read the content
            sql.or_(*[reader.matches_clause(db) for reader in readers]),

            # Conditions generated from the list of ids. These include a
            # text-ref text-content binding to connect with id data.
            *get_clauses(id_str_list, db)
            )
    else:
        readings_query = None
    return readings_query.distinct()


# =============================================================================
# Core Reading Functions
# =============================================================================


def make_db_readings(id_str_list, readers, db=None, **kwargs):
    """Read contents retrieved from the database.

    The content will be retrieved in batchs, given by the `batch` argument.
    This prevents the system RAM from being overloaded.

    Parameters
    ----------
    id_str_list : list of stings
        A list of id strings, as would retrieved from a file, each of the form
        '<id type>:<id value>', for example 'pmid:12345'
    readers : list of reader objects
        A list of the readers that will be use, for example ['reach'] if you
        wanted to use the reach reader.
    batch : int
        The number of content entries read for each batch. Default 1000.
    force_fulltext : bool
        If True, only get fulltext content from the database. Default False.

    Other keyword arguments are passed to the `read_content` function.

    Returns
    -------
    outputs : list of ReadingData instances
        The results of the readings with relevant metadata.
    """
    if db is None:
        db = get_primary_db()

    # Retrieve the kwargs needed in this function.
    batch_size = kwargs.pop('batch', 1000)
    force_fulltext = kwargs.pop('force_fulltext', False)
    force_read = kwargs.pop('force_read', False)

    # Get the iterator.
    logger.debug("Getting iterator.")
    tc_read_q = get_content_query(
        id_str_list,
        readers,
        db=db,
        force_fulltext=force_fulltext,
        force_read=force_read
        )
    logger.debug("Begginning to iterate.")
    batch_list_dict = {r.name: [] for r in readers}
    new_outputs = []
    if tc_read_q is not None:
        for text_content in tc_read_q.yield_per(batch_size):
            # The get_content function returns an iterator which yields
            # results in batches, so as not to overwhelm RAM. We need to read
            # in batches for much the same reaason.
            for r in readers:
                # Try to get a previous reading from this reader.
                reading = db.select_one(
                    db.Readings,
                    db.Readings.text_content_id == text_content.id,
                    r.matches_clause(db)
                    )
                if reading is not None and not force_read:
                    continue
                batch_list_dict[r.name].append(text_content)

                if (len(batch_list_dict[r.name])+1) % batch_size is 0:
                    # TODO: this is a bit cludgy...maybe do this better?
                    # Perhaps refactor read_content.
                    logger.debug("Reading batch of files for %s." % r.name)
                    results = r.read(batch_list_dict[r.name], **kwargs)
                    if results is not None:
                        new_outputs += results
                    batch_list_dict[r.name] = []
        logger.debug("Finished iteration.")
        # Pick up any stragglers.
        for r in readers:
            if len(batch_list_dict[r.name]) > 0:
                logger.debug("Reading remaining files for %s." % r.name)
                results = r.read(batch_list_dict[r.name], **kwargs)
                if results is not None:
                    new_outputs += results
    return new_outputs


def get_db_readings(id_lines, readers, force_fulltext=False, batch_size=1000,
                    db=None):
    """Get readings from the database."""
    if db is None:
        db = get_primary_db()

    # Get any previous readings. Note that we do this BEFORE posting the new
    # readings. Otherwise we would have duplicates.
    previous_readings_query = get_readings_query(
        id_lines,
        readers,
        db=db,
        force_fulltext=force_fulltext
        )
    if previous_readings_query is not None:
        prev_readings = [
            ReadingData(
                r.text_content_id,
                r.reader,
                r.reader_version,
                r.format,
                zlib.decompress(r.bytes, 16+zlib.MAX_WBITS).decode('utf8'),
                r.id
                )
            for r in previous_readings_query.yield_per(batch_size)
            ]
    else:
        prev_readings = []
    return prev_readings


def upload_readings(output_list, db=None):
    """Put the reading output on the database."""
    if db is None:
        db = get_primary_db()

    # Create the list of records to be copied, ensuring no uniqueness conflicts
    r_list = db.select_all(
        db.Readings,
        db.Readings.text_content_id.in_([rd.tcid for rd in output_list])
        )
    exisiting_tcid_set = set([r.text_content_id for r in r_list])
    upload_list = []
    for reading_data in output_list:
        # First check if this tcid is even in the set of existing tcids in the
        # readings table.
        if reading_data.tcid in exisiting_tcid_set:
            r_tcid_list = [r for r in r_list
                           if r.text_content_id == reading_data.tcid]
            # Now check for any exact matches:
            if any([reading_data.matches(r) for r in r_tcid_list]):
                continue

        # If there were no conflicts, we can add this to the copy list.
        upload_list.append(reading_data.make_tuple())

    # Copy into the database.
    logger.info("Adding %d/%d reading entries to the database." %
                (len(upload_list), len(output_list)))
    db.copy('readings', upload_list, ReadingData.get_cols())
    return


def produce_readings(input_list, reader_list, verbose=False, force_read=False,
                     force_fulltext=False, batch_size=1000, no_read=False,
                     no_upload=False, pickle_result=False, db=None):
    """Produce the reading output for the given ids, and upload them to db.

    This function will also retrieve pre-existing readings from the database,
    thus improving performance.

    Parameters
    ----------
    input_list : list [str]
        A list of input strings.
    reader_list : list [Reader]
        A list of Reader descendents to be used in reading.
    verbose : bool
        Optional, default False - If True, log and print the output of the
        commandline reader utilities, if False, don't.
    force_read : bool
        Optional, default False - If True, read content even if a there is an
        existing reading in the database.
    force_fulltext : bool
        Optional, default False - If True, only read fulltext article, ignoring
        abstracts.
    batch_size : int
        Optional, default 1000 - The number of text content entries to be
        yielded by the database at a given time.
    no_read : bool
        Optional, default False - If True, do not perform any new readings, and
        only retrieve existing readings from the database.
    no_upload : bool
        Optional, default False - If True, do not upload content to the
        database.
    pickle_result : bool
        Optional, default False - If True, make a pickle file of the results.
    db : indra.db.DatabaseManager instance
        Optional, default the primary database provided by `get_primary_db`
        function. Used to interface with a different databse.

    Returns
    -------
    outputs : list [ReadingData]
        A list of the outputs of the readings in the form of ReadingData
        instances.
    """
    if db is None:
        db = get_primary_db()

    prev_readings = []
    if not force_read:
        prev_readings = get_db_readings(input_list, reader_list,
                                        force_fulltext, batch_size, db=db)
    outputs = []
    if not no_read:
        outputs = make_db_readings(input_list, reader_list, verbose=verbose,
                                   force_read=force_read, db=db,
                                   force_fulltext=force_fulltext,
                                   batch=batch_size)

    if pickle_result:
        reading_out_path = pjoin(os.getcwd(), 'reading_outputs.pkl')
        with open(reading_out_path, 'wb') as f:
            pickle.dump([output.make_tuple() for output in outputs], f)
        print("Reading outputs stored in %s." % reading_out_path)

    if not no_upload:
        upload_readings(outputs, db=db)

    outputs += prev_readings

    return outputs


# =============================================================================
# Statement Processing
# =============================================================================


def upload_statements(stmt_data_list, db=None):
    """Upload the statements to the database."""
    if db is None:
        db = get_primary_db()

    logger.info("Uploading %d statements to the database." %
                len(stmt_data_list))
    db.copy('statements', [s.make_tuple() for s in stmt_data_list],
            StatementData.get_cols())

    logger.info("Uploading agents to the database.")
    reading_id_set = set([sd.reading_id for sd in stmt_data_list])
    db.insert_agents([sd.statement for sd in stmt_data_list],
                     db.Statements.reader_ref.in_(reading_id_set))
    return


def produce_statements(output_list, enrich=True, no_upload=False,
                       pickle_result=False, db=None):
    """Convert the reader output into a list of StatementData instances."""
    if db is None:
        db = get_primary_db()

    if enrich:
        _enrich_reading_data(output_list, db=db)

    stmt_data_list = make_statements(output_list)

    if not no_upload:
        upload_statements(stmt_data_list, db=db)
    if pickle_result:
        stmts_path = pjoin(os.getcwd(), 'statements.pkl')
        with open(stmts_path, 'wb') as f:
            pickle.dump([sd.statement for sd in stmt_data_list], f)
        print("Statements pickled in %s." % stmts_path)

    return stmt_data_list


# =============================================================================
# Main for script use
# =============================================================================


class ReadDBError(Exception):
    pass


if __name__ == "__main__":
    # Process the arguments. =================================================
    if args.id_file and args.file_file:
        raise ReadDBError("Cannot process both files and ids.")

    # Get the ids or files.
    with open(args.id_file, 'r') as f:
        input_lines = f.readlines()
    logger.info("Found %d ids." % len(input_lines))

    # Select only a sample of the lines, if sample is chosen.
    if args.sample is not None:
        input_lines = random.sample(input_lines, args.sample)

    # If a range is specified, only use that range.
    if args.in_range is not None:
        start_idx, end_idx = [int(n) for n in args.in_range.split(':')]
        input_lines = input_lines[start_idx:end_idx]

    # Create a single base directory
    base_dir = _get_dir('run_%s' % ('_and_'.join(args.readers)))

    # Get the readers objects.
    readers = [reader_class(base_dir=base_dir, n_proc=args.num_procs)
               for reader_class in get_readers()
               if reader_class.name.lower() in args.readers]

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    # Read everything ========================================================
    outputs = produce_readings(input_lines, readers, verbose=verbose,
                               force_read=args.force_read,
                               force_fulltext=args.force_fulltext,
                               batch_size=args.batch, no_read=args.no_read,
                               no_upload=args.no_reading_upload,
                               pickle_result=args.pickle)

    # Convert the outputs to statements ======================================
    produce_statements(outputs, no_upload=args.no_statement_upload,
                       pickle_result=args.pickle)
