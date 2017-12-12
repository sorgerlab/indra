from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.statements import *
from indra.sources.signor import SignorProcessor, _default_csv_file
from indra.db import get_primary_db
from indra.databases import hgnc_client


def by_gene_role_type(agent_id=None, agent_ns='HGNC', role=None,
                      stmt_type=None, count=1000, do_stmt_count=True, db=None):
    """Get statements from the DB by stmt type, agent, and/or agent role.

    Parameters
    ----------
    agent_id : str
        String representing the identifier of the agent from the given
        namespace. Note: if the agent namespace argument, `agent_ns`, is set
        to 'HGNC', this function will treat `agent_id` as an HGNC gene
        symbol and perform an internal lookup of the corresponding HGNC ID.
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
    do_stmt_count : bool
        Whether or not to perform an initial statement counting step to give
        more meaningful progress messages.
    db : indra.db.DatabaseManager object.
        Optionally specify a database manager that attaches to something
        besides the primary database, for example a local databse instance.

    Returns
    -------
    list of Statements from the database corresponding to the query.
    """
    if db is None:
        db = get_primary_db()

    if not (agent_id or role or stmt_type):
        raise ValueError('At least one of agent_id, role, or stmt_type '
                         'must be specified.')
    clauses = []
    if agent_id and agent_ns == 'HGNC':
        hgnc_id = hgnc_client.get_hgnc_id(agent_id)
        if not hgnc_id:
            logger.warning('Invalid gene name: %s' % agent_id)
            return []
        clauses.extend([db.Agents.db_name == 'HGNC',
                        db.Agents.db_id == hgnc_id])
    elif agent_id:
        clauses.extend([db.Agents.db_name == agent_ns,
                        db.Agents.db_id == agent_id])
    if role:
        clauses.append(db.Agents.role == role)
    if agent_id or role:
        clauses.append(db.Agents.stmt_id == db.Statements.id)
    if stmt_type:
        clauses.append(db.Statements.type == stmt_type)
    stmts = get_statements(clauses, count=count, do_stmt_count=do_stmt_count,
                           db=db)
    return stmts


def get_statements(clauses, count=1000, do_stmt_count=True, db=None):
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
    db : indra.db.DatabaseManager object.
        Optionally specify a database manager that attaches to something
        besides the primary database, for example a local database instance.

    Returns
    -------
    list of Statements from the database corresponding to the query.
    """
    if db is None:
        db = get_primary_db()

    stmts = []
    q = db.filter_query('statements', *clauses)
    if do_stmt_count:
        print("Counting statements...")
        num_stmts = q.count()
        print("Total of %d statements" % num_stmts)
    db_stmts = q.yield_per(count)
    subset = []
    total_counter = 0
    for stmt in db_stmts:
        subset.append(stmt)
        if len(subset) == count:
            stmts.extend(_stmts_from_db_list(subset))
            subset = []
        total_counter += 1
        if total_counter % count == 0:
            if do_stmt_count:
                print("%d of %d statements" % (total_counter, num_stmts))
            else:
                print("%d statements" % total_counter)

    stmts.extend(_stmts_from_db_list(subset))
    return stmts


def _stmts_from_db_list(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)

