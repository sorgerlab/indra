from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_defaults', 'get_primary_db', 'insert_agents', 'insert_pa_stmts',
           'insert_db_stmts', 'get_abstracts_by_pmids', 'get_auth_xml_pmcids',
           'get_statements_by_gene_role_type', 'get_statements',
           'make_stmts_from_db_list']

import os
import json
import logging
from os import path
from indra.databases import hgnc_client
from indra.util.get_version import get_version
from indra.util import unzip_string
from indra.statements import Complex, SelfModification, ActiveForm,\
    stmts_from_json, Conversion, Translocation
from .database_manager import DatabaseManager, IndraDatabaseError, texttypes


logger = logging.getLogger('db_util')


DEFAULTS_FILE = path.join(path.dirname(path.abspath(__file__)), 'defaults.txt')
__PRIMARY_DB = None


def get_defaults():
    "Get the default database hosts provided in the specified `DEFAULTS_FILE`."
    default_default_file = DEFAULTS_FILE
    env_key_dict = {'primary': 'INDRADBPRIMARY', 'test': 'INDRADBTEST'}
    env = os.environ
    available_keys = {k: v for k, v in env_key_dict.items() if v in env.keys()}
    if not path.exists(default_default_file) and not available_keys:
        raise IndraDatabaseError(
            "Cannot find default file or environment vars."
            )
    elif path.exists(default_default_file):
        with open(default_default_file, 'r') as f:
            defaults_raw = f.read().splitlines()
        defaults_dict = {}
        for default_line in defaults_raw:
            key, value = default_line.split('=')
            defaults_dict[key.strip()] = value.strip()
    else:
        defaults_dict = {
            purpose: env_val for purpose, my_env_key in env_key_dict.items()
            for env_key, env_val in env.items() if my_env_key == env_key
            }
    return defaults_dict


def get_primary_db(force_new=False):
    """Get a DatabaseManager instance for the primary database host.

    The primary database host is defined in the defaults.txt file, or in a file
    given by the environment variable DEFAULTS_FILE. Alternatively, it may be
    defined by the INDRADBPRIMARY environment variable. If none of the above
    are specified, this function will raise an exception.

    Note: by default, calling this function twice will return the same
    `DatabaseManager` instance. In other words:

    > db1 = get_primary_db()
    > db2 = get_primary_db()
    > db1 is db2
    True

    This means also that, for example `db1.select_one(db2.TextRef)` will work,
    in the above context.

    It is still recommended that when creating a script or function, or other
    general application, you should not rely on this feature to get your access
    to the database, as it can make substituting a different database host both
    complicated and messy. Rather, a database instance should be explicitly
    passed between different users as is done in `get_statements_by_gene_role_type`
    function's call to `get_statements` in `indra.db.query_db_stmts`.

    Parameters
    ----------
    force_new : bool
        If true, a new instance will be created and returned, regardless of
        whether there is an existing instance or not. Default is False, so that
        if this function has been called before within the global scope, a the
        instance that was first created will be returned.

    Returns
    -------
    primary_db : DatabaseManager instance
        An instance of the database manager that is attached to the primary
        database.
    """
    defaults = get_defaults()
    if 'primary' in defaults.keys():
        primary_host = defaults['primary']
    else:
        raise IndraDatabaseError("No primary host available in defaults file.")

    global __PRIMARY_DB
    if __PRIMARY_DB is None or force_new:
        __PRIMARY_DB = DatabaseManager(primary_host, label='primary')
        __PRIMARY_DB.grab_session()
    return __PRIMARY_DB


def insert_agents(db, stmt_tbl_obj, agent_tbl_obj, *other_stmt_clauses,
                  **kwargs):
    """Insert agents for statements that don't have any agents.

    Note: This method currently works for both Statements and PAStatements and their
    corresponding agents (Agents and PAAgents).

    Parameters:
    -----------
    db : indra.db.DatabaseManager
        The manager for the database into which you are adding agents.
    stmt_tbl_obj : sqlalchemy table object
        For example, `db.Statements`. The object corresponding to the
        statements column you creating agents for.
    agent_tbl_obj : sqlalchemy table object
        That agent table corresponding to the statement table above.
    *other_stmt_clauses : sqlalchemy clauses
        Further arguments, such as `db.Statements.db_ref == 1' are used to
        restrict the scope of statements whose agents may be added.
    verbose : bool
        If True, print extra information and a status bar while compiling
        agents for insert from statements. Default False.
    num_per_yield : int
        To conserve memory, statements are loaded in batches of `num_per_yeild`
        using the `yeild_per` feature of sqlalchemy queries.
    """
    verbose = kwargs.pop('verbose', False)
    num_per_yield = kwargs.pop('num_per_yield', 100)
    if len(kwargs):
        raise IndraDatabaseError("Unrecognized keyword argument(s): %s."
                                 % kwargs)
    # Build a dict mapping stmt UUIDs to statement IDs
    logger.info("Getting %s that lack %s in the database."
                % (stmt_tbl_obj.__tablename__, agent_tbl_obj.__tablename__))
    stmts_w_agents_q = db.filter_query(
        stmt_tbl_obj,
        stmt_tbl_obj.id == agent_tbl_obj.stmt_id
        )
    stmts_wo_agents_q = (db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
                         .except_(stmts_w_agents_q))
    if verbose:
        num_stmts = stmts_wo_agents_q.count()
        print("Adding agents for %d statements." % num_stmts)
    stmts_wo_agents = stmts_wo_agents_q.yield_per(num_per_yield)

    # Now assemble agent records
    logger.info("Building agent data for insert...")
    if verbose:
        print("Loading:", end='', flush=True)
    agent_data = []
    for i, db_stmt in enumerate(stmts_wo_agents):
        # Convert the database statement entry object into an indra statement.
        stmt = stmts_from_json(json.loads(db_stmt.json.decode()))

        # Figure out how the agents are structured and assign roles.
        ag_list = stmt.agent_list()
        nary_stmt_types = [Complex, SelfModification, ActiveForm, Conversion,
                           Translocation]
        if any([isinstance(stmt, tp) for tp in nary_stmt_types]):
            agents = {('OTHER', ag) for ag in ag_list}
        elif len(ag_list) == 2:
            agents = {('SUBJECT', ag_list[0]), ('OBJECT', ag_list[1])}
        else:
            raise IndraDatabaseError("Unhandled agent structure for stmt %s "
                                     "with agents: %s."
                                     % (str(stmt), str(stmt.agent_list())))

        # Prep the agents for copy into the database.
        for role, ag in agents:
            # If no agent, or no db_refs for the agent, skip the insert
            # that follows.
            if ag is None or ag.db_refs is None:
                continue
            for ns, ag_id in ag.db_refs.items():
                ag_rec = (db_stmt.id, ns, ag_id, role)
                agent_data.append(ag_rec)

        # Optionally print another tick on the progress bar.
        if verbose and i % (num_stmts//25) == 0:
            print('|', end='', flush=True)

    cols = ('stmt_id', 'db_name', 'db_id', 'role')
    db.copy(agent_tbl_obj.__tablename__, agent_data, cols)
    return


def insert_db_stmts(db, stmts, db_ref_id, verbose=False):
    """Insert statement, their database, and any affiliated agents.

    Note that this method is for uploading statements that came from a
    database to our databse, not for inserting any statements to the database.

    Parameters:
    -----------
    db : indra.db.DatabaseManager
        The manager for the database into which you are loading statements.
    stmts : list [indra.statements.Statement]
        A list of un-assembled indra statements to be uploaded to the datbase.
    db_ref_id : int
        The id to the db_ref entry corresponding to these statements.
    verbose : bool
        If True, print extra information and a status bar while compiling
        statements for insert. Default False.
    """
    # Preparing the statements for copying
    stmt_data = []
    cols = ('uuid', 'db_ref', 'type', 'json', 'indra_version')
    if verbose:
        print("Loading:", end='', flush=True)
    for i, stmt in enumerate(stmts):
        stmt_rec = (
            stmt.uuid,
            db_ref_id,
            stmt.__class__.__name__,
            json.dumps(stmt.to_json()).encode('utf8'),
            get_version()
        )
        stmt_data.append(stmt_rec)
        if verbose and i % (len(stmts)//25) == 0:
            print('|', end='', flush=True)
    if verbose:
        print(" Done loading %d statements." % len(stmts))
    db.copy('statements', stmt_data, cols)
    insert_agents(db, db.Statements, db.Agents,
                  db.Statements.db_ref == db_ref_id)
    return


def insert_pa_stmts(db, stmts, verbose=False):
    """Insert pre-assembled statements, and any affiliated agents.

    Parameters:
    -----------
    db : indra.db.DatabaseManager
        The manager for the database into which you are loading pre-assembled
        statements.
    stmts : list [indra.statements.Statement]
        A list of pre-assembled indra statements to be uploaded to the datbase.
    verbose : bool
        If True, print extra information and a status bar while compiling
        statements for insert. Default False.
    """
    logger.info("Beginning to insert pre-assembled statements.")
    stmt_data = []
    indra_version = get_version()
    cols = ('uuid', 'type', 'json', 'indra_version')
    if verbose:
        print("Loading:", end='', flush=True)
    for i, stmt in enumerate(stmts):
        stmt_rec = (
            stmt.uuid,
            stmt.__class__.__name__,
            json.dumps(stmt.to_json()).encode('utf8'),
            indra_version
        )
        stmt_data.append(stmt_rec)
        if verbose and i % (len(stmts)//25) == 0:
            print('|', end='', flush=True)
    if verbose:
        print(" Done loading %d statements." % len(stmts))
    db.copy('pa_statements', stmt_data, cols)
    insert_agents(db, db.PAStatements, db.PAAgents, verbose=verbose)
    return


def get_abstracts_by_pmids(db, pmid_list, unzip=True):
    "Get abstracts using the pmids in pmid_list."
    abst_list = db.filter_query(
        [db.TextRef, db.TextContent],
        db.TextContent.text_ref_id == db.TextRef.id,
        db.TextContent.text_type == 'abstract',
        db.TextRef.pmid.in_(pmid_list)
        ).all()
    if unzip:
        def unzip_func(s):
            return unzip_string(s.tobytes())
    else:
        def unzip_func(s):
            return s
    return [(r.pmid, unzip_func(c.content)) for (r, c) in abst_list]


def get_auth_xml_pmcids(db):
    tref_list = db.filter_query(
        [db.TextRef, db.TextContent],
        db.TextRef.id == db.TextContent.text_ref_id,
        db.TextContent.text_type == texttypes.FULLTEXT,
        db.TextContent.source == 'pmc_auth'
        )
    return [tref.pmcid for tref in tref_list]


def get_statements_by_gene_role_type(agent_id=None, agent_ns='HGNC', role=None,
                                     stmt_type=None, count=1000,
                                     do_stmt_count=True, db=None):
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
        logger.info("Counting statements...")
        num_stmts = q.count()
        logger.info("Total of %d statements" % num_stmts)
    db_stmts = q.yield_per(count)
    subset = []
    total_counter = 0
    for stmt in db_stmts:
        subset.append(stmt)
        if len(subset) == count:
            stmts.extend(make_stmts_from_db_list(subset))
            subset = []
        total_counter += 1
        if total_counter % count == 0:
            if do_stmt_count:
                logger.info("%d of %d statements" % (total_counter, num_stmts))
            else:
                logger.info("%d statements" % total_counter)

    stmts.extend(make_stmts_from_db_list(subset))
    return stmts


def make_stmts_from_db_list(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)
