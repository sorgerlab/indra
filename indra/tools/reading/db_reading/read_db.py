"""This module provides essential tools to run reading using indra's own
database. This may also be run as a script; for details run:
`python read_pmids_db --help`
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import random
import pickle
from math import log10, floor, ceil
from datetime import datetime
from indra.tools.reading.util.script_tools import get_parser, make_statements,\
                                             StatementData

logger = logging.getLogger('make_db_readings')
if __name__ == '__main__':
    parser = get_parser(
        'A tool to read and process content from the database.',
        ('A file containing a list of ids of the form <id_type>:<id>. '
         'Note that besided the obvious id types (pmid, pmcid, doi, etc.), '
         'you may use trid and tcid to indicate text ref and text content '
         'ids, respectively. Note that these are specific to the database, '
         'and should thus be used with care.')
        )
    parser.add_argument(
        '-m', '--reading_mode',
        choices=['all', 'unread', 'none'],
        default='unread',
        help=("Set the reading mode. If 'all', read everything, if "
              "'unread', only read content that does not have pre-existing "
              "readings of the same reader and version, if 'none', only "
              "use pre-existing readings. Default is 'unread'.")
        )
    parser.add_argument(
        '-s', '--stmt_mode',
        choices=['all', 'unread', 'none'],
        default='all',
        help=("Choose which readings should produce statements. If 'all', all "
              "readings that are produced or retrieved will be used to produce "
              "statements. If 'unread', only produce statements from "
              "previously unread content. If 'none', do not produce any "
              "statements (only readings will be produced).")
        )
    parser.add_argument(
        '-t', '--temp',
        default='.',
        help='Select the location of the temp file.'
        )
    parser.add_argument(
        '-o', '--output',
        dest='name',
        help=('Pickle all results and save in files labelled as '
              '<NAME>_<output_type>.pkl.'),
        default=None
        )
    parser.add_argument(
        '-b', '--inner_batch',
        dest='b_in',
        help=('Choose the size of the inner batches, which is the number of '
              'text content entires loaded at a given time, and the number of '
              'entries that are read at a time by a reader. The default is '
              '1,000.'),
        default=1000,
        type=int
        )
    parser.add_argument(
        '-B', '--outer_batch',
        dest='b_out',
        default=10000,
        type=int,
        help=('Select the number of ids to read per outer level batch. This '
              'determines the number of readings/statements uploaded/pickled '
              'at a time, and thus also limits the amount of RAM that will be '
              'used. A larger outer batch means more RAM. The default is '
              '10,000.')
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
        '--force_fulltext',
        help='Make the reader only read full text from the database.',
        action='store_true'
        )
    parser.add_argument(
        '--use_best_fulltext',
        help='Use only the best full text available.',
        action='store_true'
        )
    parser.add_argumet(
        '--max_reach_space_ratio',
        type=float,
        help='Set the maximum ratio of spaces to non-spaces for REACH input.',
        default=None
    )
    parser.add_argument(
        '--max_reach_input_len',
        type=int,
        help='Set the maximum length of content that REACH will read.',
        default=None
    )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

from indra.literature.elsevier_client import extract_text as process_elsevier
from indra.db import get_primary_db, formats, texttypes
from indra.db import sql_expressions as sql
from indra.db.util import insert_agents
from indra.tools.reading.readers import ReadingData, _get_dir, get_reader, \
    Content


class ReadDBError(Exception):
    pass


# =============================================================================
# Useful functions
# =============================================================================


def _convert_id_entry(id_entry, allowed_types=None):
    ret = [s.strip() for s in id_entry.split(':')]
    if len(ret) != 2:
        raise ReadDBError('Improperly formatted id entry: \"%s\"' % id_entry)
    ret[0] = ret[0].lower()
    if allowed_types is not None and ret[0] not in allowed_types:
        raise ReadDBError('Invalid id type: \"%s\"' % ret[0])
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
        'reading',
        db.Reading.text_content_id.in_([rd.tcid for rd in reading_data_iter
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


def get_id_dict(id_str_list):
    """Parse the list of id string into a dict."""
    id_types = get_primary_db().TextRef.__table__.columns.keys()
    id_types.remove('id')
    id_types += ['trid', 'tcid']
    id_dict = {id_type: [] for id_type in id_types}
    for id_entry in id_str_list:
        id_type, id_val = _convert_id_entry(id_entry, id_types)
        if id_type in ['trid', 'tcid']:
            id_dict[id_type].append(int(id_val))
        else:
            id_dict[id_type].append(id_val)
    return id_dict


def get_priority_tcids(id_dict, priorities, always_add=None, db=None):
    """For all ids, besides tcids, choose best content available.

    This function will convert all ids to tcids.
    """
    if db is None:
        db = get_primary_db()

    def is_better(new, old):
        if new in priorities and old in priorities:
            return priorities.index(new) < priorities.index(old)
        return False

    logger.debug("Getting content prioritized by %s." % str(priorities))
    tcids = set(id_dict.pop('tcid', []))
    clauses = get_clauses(id_dict, db)
    tcid_source = set()
    for clause in clauses:
        q = (db.session.query(db.TextRef.id, db.TextContent.id,
                              db.TextContent.source)
             .filter(db.TextContent.text_ref_id == db.TextRef.id, clause))
        id_set = set(q.all())
        logger.debug("Got %d more ids." % len(id_set))
        tcid_source |= id_set
    logger.debug("Got %d id's total." % len(tcid_source))
    tr_best = {}
    for trid, tcid, source in tcid_source:
        if trid not in tr_best.keys() or is_better(source, tr_best[trid][0]):
            tr_best[trid] = (source, tcid)
        if always_add is not None and source in always_add:
            tcids.add(tcid)
    tcids |= {tcid for _, tcid in tr_best.values()}
    return tcids


def get_clauses(id_dict, db):
    """Get a list of clauses to be passed to a db query.

    Note that an empty condition will be returned if id_dict has no ids in it
    (either the dict is empty or all the lists within the dict are empty),
    which will in general have the unexpected effect of selecting everything,
    rather than nothing.

    Parameters
    ----------
    id_dict : dict {id_type: [int or str]}
        A dictionary indexed by the type of id, containing lists of id's of
        that the respective type. If all the lists are empty, or the dict is
        empty, returns an empty condition. Note that id types of 'trid' and
        'tcid' will be mapped to text ref ids and text content ids,
        respectively.
    db : indra.db.DatabaseManager instance
        This instance is only used for forming the query, and will not be
        accessed or queried.

    Returns
    -------
    clause_list : list [sqlalchemy clauses]
        A list of sqlalchemy clauses to be used in query in the form:
        `db.filter_query(<table>, <other clauses>, *clause_list)`.
        If the id_dict has no ids, an effectively empty condition is returned.
    """
    # Handle all id types but text ref ids (trid) and text content ids (tcid).
    id_condition_list = [getattr(db.TextRef, id_type).in_(id_list)
                         for id_type, id_list in id_dict.items()
                         if len(id_list) and id_type not in ['tcid', 'trid']]

    # Handle the special id types trid and tcid.
    for id_type, table in [('trid', db.TextRef), ('tcid', db.TextContent)]:
        if id_type in id_dict.keys() and len(id_dict[id_type]):
            int_id_list = [int(i) for i in id_dict[id_type]]
            id_condition_list.append(table.id.in_(int_id_list))
    return id_condition_list


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


def get_content_query(ids, readers, db=None, force_fulltext=False,
                      force_read=False, debug=False, print_summary=False):
    """Construct a query to access all the content that will be read.

    If ids is not 'all', and does not contain any ids, None is returned.

    Parameters
    ----------
    ids : 'all' or dict {<id type> : [str/int]}
        If 'all', then all the content will be included in the query. Otherwise
        a the content will be constrained to that corresponding to the ids in
        id_dict, which are matched using text refs.
    readers : list [Reader child instances]
        A list of the reader objects, which contain the required metadata (name
        and version of the reader) used to find content that needs to be read.
    db : indra.db.DatabaseManager instance
        Optional, default None, in which case the primary database is used. If
        specified, the alternative database will be used. This function should
        not alter the database.
    force_fulltext : bool
        Optional, default False - If True, only fulltext content will be read,
        as opposed to including abstracts.
    force_read : bool
        Optional, default False - If True, all content will be returned,
        whether it has been read or not.

    Returns
    -------
    tc_tbr_query : sqlalchemy query object or None
        The query of the text content to be read (tc_tbr). If there are no ids
        contained in ids, or it is not 'all', return None.
    """
    if debug:
        logger.setLevel(logging.DEBUG)
    if db is None:
        db = get_primary_db()
    logger.debug("Got db handle.")

    # These allow conditions on different tables to equal conditions on the
    # dependent tables.
    tc_tr_binding = db.TextContent.text_ref_id == db.TextRef.id
    rd_tc_binding = db.Reading.text_content_id == db.TextContent.id

    # Begin the list of clauses with the binding between text content and
    # text refs.
    clauses = [tc_tr_binding]

    # Add a fulltext requirement, if applicable.
    if force_fulltext:
        clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)

    # If we are actually getting anything, else we return None.
    if ids == 'all' or any([len(id_list) > 0 for id_list in ids.values()]):
        if ids is not 'all':
            sub_clauses = get_clauses(ids, db)
            if len(sub_clauses) > 1:
                clauses.append(sql.or_(*sub_clauses))
            else:
                clauses.append(*sub_clauses)

        # Get the text content query object
        tc_query = db.filter_query(
            db.TextContent,
            *clauses
            ).distinct()

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

        if print_summary:
            try:
                logger.debug("Going to try to make a nice summary...")
                logger.info(get_text_content_summary_string(tc_tbr_query, db))
            except Exception:
                logger.debug("Could not print summary of results.")
    else:
        logger.debug("No ids in id_dict, so no query formed.")
        return None

    return tc_tbr_query.distinct()


def get_readings_query(ids, readers, db=None, force_fulltext=False):
    """Create a query to access all the relevant existing readings.

    Note that if ids is not 'all' and ids is a dict with no ids in it,
    this function returns None.

    Parameters
    ----------
    ids : 'all' or dict {<id_type> : [str/int]}
        If 'all', then all possible readings in the database matching the given
        readers and other conditions will be returned. Otherwise, only those
        that correspond to one of the ids in ids dict will be contained. If an
        ids dict has no ids in it, None is returned.
    readers : list [Reader child instances]
        A list of the readers whose names and versions you wish to match in the
        readings queried from the database.
    db : indra.db.DatabaseManager instance
        Optional, default None, in which case the primary database is used. If
        specified, the alternative database will be used. This function should
        not alter the database.
    force_fulltext : bool
        Optional, default False - If True, only readings corresponding to
        fulltext content will be read, as opposed to including readings created
        from abstracts.

    Returns
    -------
    readings_query : sql query instance or None
        Returns a query that can be used to access the specified content, or
        else None if no content was specified.
    """
    if db is None:
        db = get_primary_db()
    clauses = [
        # Bind conditions on readings to conditions on content.
        db.Reading.text_content_id == db.TextContent.id,

        # Bind text content to text refs
        db.TextContent.text_ref_id == db.TextRef.id,

        # Check if at least one of the readers has read the content
        sql.or_(*[reader.matches_clause(db) for reader in readers])
        ]
    if force_fulltext:
        clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)

    if ids == 'all' or any([id_list for id_list in ids.values()]):
        if ids != 'all':
            sub_clauses = get_clauses(ids, db)
            if len(sub_clauses) > 1:
                clauses.append(sql.or_(*sub_clauses))
            else:
                clauses.append(*sub_clauses)

        readings_query = db.filter_query(
            db.Reading,

            # Bind conditions on readings to conditions on content.
            db.Reading.text_content_id == db.TextContent.id,

            # Bind text content to text refs
            db.TextContent.text_ref_id == db.TextRef.id,

            # Check if at least one of the readers has read the content
            sql.or_(*[reader.matches_clause(db) for reader in readers]),

            # Conditions generated from the list of ids. These include a
            # text-ref text-content binding to connect with id data.
            *clauses
            )
    else:
        return None

    return readings_query.distinct()


def process_content(text_content):
    """Get the appropriate content object from the text content."""
    if text_content.format == formats.TEXT:
        cont_fmt = 'txt'
    elif (text_content.source in ['pmc_oa', 'manuscripts']
          and text_content.format == formats.XML):
        cont_fmt = 'nxml'
    else:
        cont_fmt = text_content.format
    content = Content.from_string(text_content.id, cont_fmt,
                                  text_content.content, compressed=True,
                                  encoded=True)
    if text_content.source == 'elsevier':
        raw_xml_text = content.get_text()
        elsevier_text = process_elsevier(raw_xml_text)
        if elsevier_text is None:
            logger.warning("Could not extract text from Elsevier xml for tcid: "
                           "%d" % text_content.id)
            return None
        content = Content.from_string(content.get_id(), 'text', elsevier_text)
    return content


# =============================================================================
# Core Reading Functions
# =============================================================================


def make_db_readings(id_dict, readers, batch_size=1000, force_fulltext=False,
                     force_read=False, skip_dict=None, db=None, **kwargs):
    """Read contents retrieved from the database.

    The content will be retrieved in batchs, given by the `batch` argument.
    This prevents the system RAM from being overloaded.

    Parameters
    ----------
    id_dict : dict {<id_type>:[<id value>, ...]}
        A dict of lists of the id's to be read, keyed by id_type.
    readers : list of reader objects
        A list of the readers that will be use, for example ['reach'] if you
        wanted to use the reach reader.
    batch_size : int
        The number of content entries read for each batch. Default 1000.
    force_fulltext : bool
        If True, only get fulltext content from the database. Default False.
    force_read : bool
        If True, read even if text_content id is found in skip_dict.
    skip_dict : dict {<reader> : list [int]}
        A dict containing text content id's to be skipped.
    db : indra.db.DatabaseManager instance
        A handle to a database. Default None; if None, a handle to the primary
        database (see indra.db) is retrieved.

    Other keyword arguments are passed to the `read` methods of the readers.

    Returns
    -------
    outputs : list of ReadingData instances
        The results of the readings with relevant metadata.
    """
    if db is None:
        db = get_primary_db()

    # Get the iterator.
    logger.debug("Getting iterator.")
    tc_read_q = get_content_query(
        id_dict,
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
            # in batches for much the same reason.
            for r in readers:
                if not force_read:
                    if skip_dict is not None:
                        if text_content.id in skip_dict[r.name]:
                            continue
                    else:
                        # Try to get a previous reading from this reader.
                        reading = db.select_one(
                            db.Reading,
                            db.Reading.text_content_id == text_content.id,
                            r.matches_clause(db)
                            )
                        if reading is not None:
                            continue
                processed_content = process_content(text_content)
                if processed_content is not None:
                    batch_list_dict[r.name].append(processed_content)

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


def get_db_readings(id_dict, readers, force_fulltext=False, batch_size=1000,
                    db=None):
    """Get readings from the database."""
    if db is None:
        db = get_primary_db()

    # Get any previous readings. Note that we do this BEFORE posting the new
    # readings. Otherwise we would have duplicates.
    previous_readings_query = get_readings_query(
        id_dict,
        readers,
        db=db,
        force_fulltext=force_fulltext
        )
    if previous_readings_query is not None:
        prev_readings = [
            ReadingData.from_db_reading(r)
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
        db.Reading,
        db.Reading.text_content_id.in_([rd.tcid for rd in output_list])
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
    db.copy('reading', upload_list, ReadingData.get_cols())
    return


def produce_readings(id_dict, reader_list, verbose=False, read_mode='unread',
                     get_preexisting=True, force_fulltext=False,
                     batch_size=1000, no_upload=False, pickle_file=None,
                     db=None, log_readers=True, prioritize=False):
    """Produce the reading output for the given ids, and upload them to db.

    This function will also retrieve pre-existing readings from the database,
    thus improving performance.

    Parameters
    ----------
    id_dict : dict {<id_type>:[<id value>, ...]}
        A dict of lists of the id's to be read, keyed by id_type.
    reader_list : list [Reader]
        A list of Reader descendents to be used in reading.
    verbose : bool
        Optional, default False - If True, log and print the output of the
        commandline reader utilities, if False, don't.
    read_mode : str : 'all', 'unread', or 'none'
        Optional, default 'undread' - If 'all', read everything (generally
        slow); if 'unread', only read things that were unread, (the cache of old
        readings may still be used if `stmt_mode='all'` to get everything); if
        'none', don't read, and only retrieve existing readings.
    get_preexisting : bool
        Optional, default True. If True, retrieve old readings where available
        (if `read_mode` is not 'all'). If False, don't retrieve old readings.
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
    pickle_file : str or None
        Optional, default None - otherwise the path to a file in which the
        reading data will be saved.
    db : indra.db.DatabaseManager instance
        Optional, default is None, in which case the primary database provided
        by `get_primary_db` function is used. Used to interface with a
        different databse.
    log_readers : bool
        Default True. If True, stash the logs of the readers in a file.
    prioritize : bool
        Default False. If True, choose only the best content to read.

    Returns
    -------
    outputs : list [ReadingData]
        A list of the outputs of the readings in the form of ReadingData
        instances.
    """
    # Get a database instance.
    logger.debug("Producing readings in %s mode." % read_mode)
    if db is None:
        db = get_primary_db()

    # Sort out our priorities
    if prioritize:
        logger.debug("Prioritizing...")
        tcids = get_priority_tcids(id_dict,
                                   ['pmc_oa', 'manuscripts', 'elsevier'],
                                   always_add=['pubmed'], db=db)
        id_dict = {'tcid': list(tcids)}

    # Handle the cases where I need to retrieve old readings.
    prev_readings = []
    skip_reader_tcid_dict = None
    if get_preexisting and read_mode != 'all':
        prev_readings = get_db_readings(id_dict, reader_list, force_fulltext,
                                        batch_size, db=db)
        skip_reader_tcid_dict = {r.name: [] for r in reader_list}
        logger.info("Found %d pre-existing readings." % len(prev_readings))
        if read_mode != 'none':
            for rd in prev_readings:
                skip_reader_tcid_dict[rd.reader].append(rd.tcid)

    # Now produce any new readings that need to be produced.
    outputs = []
    if read_mode != 'none':
        outputs = make_db_readings(id_dict, reader_list, verbose=verbose,
                                   skip_dict=skip_reader_tcid_dict, db=db,
                                   force_fulltext=force_fulltext,
                                   force_read=(read_mode == 'all'),
                                   batch_size=batch_size, log=log_readers)
        logger.info("Made %d new readings." % len(outputs))

    if not no_upload:
        try:
            upload_readings(outputs, db=db)
        except Exception as e:
            logger.exception(e)
            if pickle_file is None:
                pickle_file = ("failure_reading_dump_%s.pkl"
                               % datetime.now().strftime('%Y%m%d_%H%M%S'))
            logger.error("Cound not upload readings. Results are pickled in: "
                         + pickle_file)

    outputs += prev_readings

    if pickle_file is not None:
        with open(pickle_file, 'wb') as f:
            pickle.dump([output.make_tuple() for output in outputs], f)
        print("Reading outputs stored in %s." % pickle_file)

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
    db.copy('raw_statements', [s.make_tuple() for s in stmt_data_list],
            StatementData.get_cols())

    logger.info("Uploading agents to the database.")
    reading_id_set = set([sd.reading_id for sd in stmt_data_list])
    if len(reading_id_set):
        db_stmts = (
            db.select_one(db.RawStatements,
                          db.RawStatements.uuid.like(s.statement.uuid))
            for s in stmt_data_list
            )
        insert_agents(db, 'raw', db_stmts, verbose=True)
    return


def produce_statements(output_list, enrich=True, no_upload=False,
                       pickle_file=None, n_proc=1, db=None):
    """Convert the reader output into a list of StatementData instances."""
    if db is None:
        db = get_primary_db()

    if enrich:
        _enrich_reading_data(output_list, db=db)

    stmt_data_list = make_statements(output_list, n_proc)

    if not no_upload:
        try:
            upload_statements(stmt_data_list, db=db)
        except Exception as e:
            logger.exception(e)
            if pickle_file is None:
                pickle_file = ("failure_stmt_dump_%s.pkl"
                               % datetime.now().strftime('%Y%m%d_%H%M%S'))
            logger.error("Could not upload statements. Results pickled in: %s."
                         % pickle_file)
    if pickle_file is not None:
        with open(pickle_file, 'wb') as f:
            pickle.dump([sd.statement for sd in stmt_data_list], f)
        print("Statements pickled in %s." % pickle_file)

    return stmt_data_list


# =============================================================================
# Main for script use
# =============================================================================


if __name__ == "__main__":
    # Process the arguments. =================================================

    # Get the ids.
    with open(args.input_file, 'r') as f:
        input_lines = f.readlines()
    logger.info("Found %d ids." % len(input_lines))

    # Select only a sample of the lines, if sample is chosen.
    if args.n_samp is not None:
        input_lines = random.sample(input_lines, args.n_samp)
    else:
        random.shuffle(input_lines)

    # If a range is specified, only use that range.
    if args.range_str is not None:
        start_idx, end_idx = [int(n) for n in args.range_str.split(':')]
        input_lines = input_lines[start_idx:end_idx]

    # Get the outer batch.
    B = args.b_out
    n_max = int(ceil(float(len(input_lines))/B))

    # Create a single base directory
    base_dir = _get_dir(args.temp, 'run_%s' % ('_and_'.join(args.readers)))

    # Get the readers objects.
    special_reach_args_dict = {
        'input_character_limit': args.max_reach_space_ratio,
        'max_space_ratio': args.max_reach_input_len
    }
    readers = []
    for reader_name in args.readers:
        kwargs = {'base_dir': base_dir, 'n_proc': args.n_proc}
        if reader_name == 'REACH':
            for key_name, reach_arg in special_reach_args_dict.items():
                if reach_arg is not None:
                    kwargs[key_name] = reach_arg
        readers.append(get_reader(reader_name, **kwargs))

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    # Some combinations of options don't make sense:
    forbidden_combos = [('all', 'unread'), ('none', 'unread'), ('none', 'none')]
    assert (args.reading_mode, args.stmt_mode) not in forbidden_combos, \
        ("The combination of reading mode %s and statement mode %s is not "
         "allowed." % (args.reading_mode, args.stmt_mode))

    for n in range(n_max):
        logger.info("Beginning outer batch %d/%d. ------------" % (n+1, n_max))

        # Get the pickle file names.
        if args.name is not None:
            reading_pickle = args.name + '_readings_%d.pkl' % n
            stmts_pickle = args.name + '_stmts_%d.pkl' % n
        else:
            reading_pickle = None
            stmts_pickle = None

        # Get the dict of ids.
        id_dict = get_id_dict(input_lines[B*n:B*(n+1)])

        # Read everything ====================================================
        outputs = produce_readings(id_dict, readers, verbose=verbose,
                                   read_mode=args.reading_mode,
                                   get_preexisting=(args.stmt_mode == 'all'),
                                   batch_size=args.b_in,
                                   force_fulltext=args.force_fulltext,
                                   no_upload=args.no_reading_upload,
                                   pickle_file=reading_pickle,
                                   prioritize=args.use_best_fulltext)

        # Convert the outputs to statements ==================================
        if args.stmt_mode != 'none':
            produce_statements(outputs, no_upload=args.no_statement_upload,
                               pickle_file=stmts_pickle, n_proc=args.n_proc)
