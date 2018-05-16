from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
from collections import defaultdict

from indra.databases import hgnc_client
from indra.db import util as db_util

logger = logging.getLogger('db_client')


def get_reader_output(db, ref_id, ref_type='tcid', reader=None,
                      reader_version=None):
    """Return reader output for a given text content.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        Reference to the DB to query
    ref_id : int or str
        The text reference ID whose reader output should be returned
    ref_type : Optional[str]
        The type of ID to look for, options include
        'tcid' for the database's internal unique text content ID,
        or 'pmid', 'pmcid', 'doi, 'pii', 'manuscript_id'
        Default: 'tcid'
    reader : Optional[str]
        The name of the reader whose output is of interest
    reader_version : Optional[str]
        The specific version of the reader

    Returns
    -------
    reading_results : dict{dict{list[str]}}
        A dict of reader outputs that match the query criteria, indexed first
        by text content id, then by reader.
    """
    if ref_type == 'tcid':
        clauses = [db.Readings.text_content_id == ref_id]
    else:
        trids = _get_trids(db, ref_id, ref_type)
        if not trids:
            return []
        logger.debug("Found %d text ref ids." % len(trids))
        clauses = [db.TextContent.text_ref_id.in_(trids),
                   db.Readings.text_content_id == db.TextContent.id]
    if reader:
        clauses.append(db.Readings.reader == reader.upper())
    if reader_version:
        clauses.append(db.Readings.reader_version == reader_version)

    res = db.select_all([db.Readings.text_content_id, db.Readings.reader,
                         db.Readings.bytes], *clauses)
    reading_dict = defaultdict(lambda: defaultdict(lambda: []))
    for tcid, reader, result in res:
        reading_dict[tcid][reader].append(db_util.unpack(result))
    return reading_dict


def get_content_by_refs(db, pmid_list=None, trid_list=None, sources=None,
                        formats=None, content_type='abstract', unzip=True):
    """Return content from the database given a list of PMIDs or text ref ids.

    Note that either pmid_list OR trid_list must be set, and only one can be
    set at a time.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        Reference to the DB to query
    pmid_list : list[str] or None
        A list of pmids. Default is None, in which case trid_list must be given.
    trid_list : list[int] or None
        A list of text ref ids. Default is None, in which case pmid list must be
        given.
    sources : list[str] or None
        A list of sources to include (e.g. 'pmc_oa', or 'pubmed'). Default is
        None, indicating that all sources will be included.
    formats : list[str]
        A list of the formats to be included ('xml', 'text'). Default is None,
        indicating that all formats will be included.
    content_type : str
        Select the type of content to load ('abstract' or 'fulltext'). Note that
        not all refs will have any, or both, types of content.
    unzip : Optional[bool]
        If True, the compressed output is decompressed into clear text.
        Default: True

    Returns
    -------
    content_dict : dict
        A dictionary whose keys are text ref ids, with each value being the
        the corresponding content.
    """
    # Make sure we only get one type of list.
    if not pmid_list or trid_list:
        raise ValueError("One of `pmid_list` or `trid_list` must be defined.")
    if pmid_list and trid_list:
        raise ValueError("Only one of `pmid_list` or `trid_list` may be used.")

    # Put together the clauses for the general constraints.
    clauses = []
    if sources is not None:
        clauses.append(db.TextContent.source.in_(sources))
    if formats is not None:
        clauses.append(db.TextContent.format.in_(formats))
    if content_type not in ['abstract', 'fulltext']:
        raise ValueError("Unrecognized content type: %s" % content_type)
    else:
        clauses.append(db.TextContent.text_type == content_type)

    # Do the query to get the content.
    if pmid_list is not None:
        content_list = db.select_all([db.TextRef.pmid, db.TextContent.content],
                                     db.TextRef.id == db.TextContent.text_ref_id,
                                     db.TextRef.pmid.in_(pmid_list),
                                     *clauses)
    else:
        content_list = db.select_all([db.TextRef.id, db.TextContent.content],
                                     db.TextContent.text_ref_id.in_(trid_list),
                                     *clauses)
    if unzip:
        content_dict = {id_val: db_util.unpack(content)
                        for id_val, content in content_list}
    else:
        content_dict = {id_val: content for id_val, content in content_list}
    return content_dict


#==============================================================================
# Below are some functions that are useful for getting raw statements from the
# database at various levels of abstraction.
#==============================================================================

def get_statements_by_gene_role_type(agent_id=None, agent_ns='HGNC-SYMBOL',
                                     role=None, stmt_type=None, count=1000,
                                     db=None, do_stmt_count=True,
                                     preassembled=True):
    """Get statements from the DB by stmt type, agent, and/or agent role.

    Parameters
    ----------
    agent_id : str
        String representing the identifier of the agent from the given
        namespace. Note: if the agent namespace argument, `agent_ns`, is set
        to 'HGNC-SYMBOL', this function will treat `agent_id` as an HGNC gene
        symbol and perform an internal lookup of the corresponding HGNC ID.
        Default is 'HGNC-SYMBOL'.
    agent_ns : str
        Namespace for the identifier given in `agent_id`.
    role : str
        String corresponding to the role of the agent in the statement.
        Options are 'SUBJECT', 'OBJECT', or 'OTHER' (in the case of `Complex`,
        `SelfModification`, and `ActiveForm` Statements).
    stmt_type : str
        Name of the Statement class.
    count : int
        Number of statements to retrieve in each batch (passed to
        :py:func:`get_statements`).
    db : :py:class:`DatabaseManager`
        Optionally specify a database manager that attaches to something
        besides the primary database, for example a local databse instance.
    do_stmt_count : bool
        Whether or not to perform an initial statement counting step to give
        more meaningful progress messages.
    preassembled : bool
        If true, statements will be selected from the table of pre-assembled
        statements. Otherwise, they will be selected from the raw statements.
        Default is True.

    Returns
    -------
    list of Statements from the database corresponding to the query.
    """
    if db is None:
        db = db_util.get_primary_db()

    if preassembled:
        Statements = db.PAStatements
        Agents = db.PAAgents
    else:
        Statements = db.Statements
        Agents = db.Agents

    if not (agent_id or role or stmt_type):
        raise ValueError('At least one of agent_id, role, or stmt_type '
                         'must be specified.')
    clauses = []
    if agent_id and agent_ns == 'HGNC-SYMBOL':
        hgnc_id = hgnc_client.get_hgnc_id(agent_id)
        if not hgnc_id:
            logger.warning('Invalid gene name: %s' % agent_id)
            return []
        clauses.extend([Agents.db_name.like('HGNC'),
                        Agents.db_id.like(hgnc_id)])
    elif agent_id:
        clauses.extend([Agents.db_name.like(agent_ns),
                        Agents.db_id.like(agent_id)])
    if role:
        clauses.append(Agents.role == role)
    if agent_id or role:
        clauses.append(Agents.stmt_id == Statements.id)
    if stmt_type:
        clauses.append(Statements.type == stmt_type)
    stmts = get_statements(clauses, count=count, do_stmt_count=do_stmt_count,
                           db=db, preassembled=preassembled)
    return stmts


def get_statements_by_paper(id_val, id_type='pmid', count=1000, db=None,
                            do_stmt_count=True):
    """Get the statements from a particular paper.

    Note: currently this can only retrieve raw statements, because of the
    partially implemented configuration of the pre-assembled Statement table.

    Parameters
    ----------
    id_val : int or str
        The value of the id for the paper whose statements you wish to retrieve.
    id_type : str
        The type of id used (default is pmid). Options include pmid, pmcid, doi,
        pii, url, or manuscript_id. Note that pmid is generally the best means
        of getting a paper.
    count : int
        Number of statements to retrieve in each batch (passed to
        :py:func:`get_statements`).
    db : :py:class:`DatabaseManager`
        Optionally specify a database manager that attaches to something
        besides the primary database, for example a local databse instance.
    do_stmt_count : bool
        Whether or not to perform an initial statement counting step to give
        more meaningful progress messages.

    Returns
    -------
    A list of Statements from the database corresponding to the paper id given.
    """
    if db is None:
        db = db_util.get_primary_db()

    trid_list = _get_trids(db, id_val, id_type)
    if not trid_list:
        return None

    stmts = []
    for trid in trid_list:
        clauses = [
            db.TextContent.text_ref_id == trid,
            db.Readings.text_content_id == db.TextContent.id,
            db.Statements.reader_ref == db.Readings.id
        ]
        stmts.extend(get_statements(clauses, count=count, preassembled=False,
                                    do_stmt_count=do_stmt_count, db=db))
    return stmts


def get_statements(clauses, count=1000, do_stmt_count=True, db=None,
                   preassembled=True):
    """Select statements according to a given set of clauses.

    Parameters
    ----------
    clauses : list
        list of sqlalchemy WHERE clauses to pass to the filter query.
    count : int
        Number of statements to retrieve and process in each batch.
    do_stmt_count : bool
        Whether or not to perform an initial statement counting step to give
        more meaningful progress messages.
    db : :py:class:`DatabaseManager`
        Optionally specify a database manager that attaches to something
        besides the primary database, for example a local database instance.
    preassembled : bool
        If true, statements will be selected from the table of pre-assembled
        statements. Otherwise, they will be selected from the raw statements.
        Default is True.

    Returns
    -------
    list of Statements from the database corresponding to the query.
    """
    if db is None:
        db = db_util.get_primary_db()

    stmts_tblname = 'pa_statements' if preassembled else 'statements'

    stmts = []
    q = db.filter_query(stmts_tblname, *clauses)
    if do_stmt_count:
        logger.info("Counting statements...")
        num_stmts = q.count()
        logger.info("Total of %d statements" % num_stmts)
    db_stmts = q.yield_per(count)
    subset = []
    total_counter = 0
    for stmt in db_stmts:
        subset.append(stmt)
        if len(subset) == count:
            stmts.extend(db_util.make_stmts_from_db_list(subset))
            subset = []
        total_counter += 1
        if total_counter % count == 0:
            if do_stmt_count:
                logger.info("%d of %d statements" % (total_counter, num_stmts))
            else:
                logger.info("%d statements" % total_counter)

    stmts.extend(db_util.make_stmts_from_db_list(subset))
    return stmts


def _get_trids(db, id_val, id_type):
    """Return text ref IDs corresponding to any ID type and value."""
    # Get the text ref id(s)
    if id_type in ['trid']:
        trid_list = [int(id_val)]
    else:
        id_types = ['pmid', 'pmcid', 'doi', 'pii', 'url', 'manuscript_id']
        if id_type not in id_types:
            raise ValueError('id_type must be one of: %s' % str(id_types))
        constraint = (getattr(db.TextRef, id_type) == id_val)
        trid_list = [trid for trid, in db.select_all(db.TextRef.id, constraint)]
    return trid_list
