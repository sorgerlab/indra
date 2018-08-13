from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from datetime import datetime

__all__ = ['get_defaults', 'get_primary_db', 'get_db', 'insert_agents',
           'insert_pa_stmts', 'insert_db_stmts', 'get_raw_stmts_frm_db_list',
           'distill_stmts']

import re
import json
import zlib
import logging
from itertools import groupby
from functools import partial, wraps
from multiprocessing.pool import Pool

from indra.util import batch_iter
from indra.util.nested_dict import NestedDict
from indra.util.get_version import get_version
from indra.statements import Complex, SelfModification, ActiveForm,\
    stmts_from_json, Conversion, Translocation, Statement

from .database_manager import DatabaseManager, IndraDatabaseError
from .config import get_databases as get_defaults

logger = logging.getLogger('db_util')


__PRIMARY_DB = None


def _clockit(func):
    @wraps(func)
    def timed_func(*args, **kwargs):
        start = datetime.now()
        ret = func(*args, **kwargs)
        end = datetime.now()
        # print(u'\033[0;35;40m%s \033[1;36;40m%-30s\033[0;35;40m %s %s \033[0m'
        #       % ('~'*5, func.__name__, end-start, '~'*5))
        print('%s %-30s %s %s' % ('~'*5, func.__name__, end-start, '~'*5))
        #fname = '%s-%s_times.log' % (abspath(__file__), func.__name__)
        #with open(fname, 'a') as f:
        #    f.write('%s: %s\n' % (start, end-start))
        return ret
    return timed_func


def get_primary_db(force_new=False):
    """Get a DatabaseManager instance for the primary database host.

    The primary database host is defined in the defaults.txt file, or in a file
    given by the environment variable DEFAULTS_FILE. Alternatively, it may be
    defined by the INDRADBPRIMARY environment variable. If none of the above
    are specified, this function will raise an exception.

    Note: by default, calling this function twice will return the same
    `DatabaseManager` instance. In other words::

        db1 = get_primary_db()
        db2 = get_primary_db()
        db1 is db2

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
    primary_db : :py:class:`DatabaseManager`
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


def get_test_db():
    """Get a DatabaseManager for the test database."""
    defaults = get_defaults()
    test_defaults = {k: v for k, v in defaults.items() if 'test' in k}
    key_list = list(test_defaults.keys())
    key_list.sort()
    db = None
    for k in key_list:
        test_name = test_defaults[k]
        m = re.match('(\w+)://.*?/([\w.]+)', test_name)
        if m is None:
            logger.warning("Poorly formed db name: %s" % test_name)
            continue
        sqltype = m.groups()[0]
        try:
            db = DatabaseManager(test_name, sqltype=sqltype, label=k)
            db.grab_session()
        except Exception as e:
            logger.error("%s didn't work" % test_name)
            logger.exception(e)
            continue  # Clearly this test database won't work.
        logger.info("Using test database %s." % k)
        break
    if db is None:
        logger.error("Could not find any test database names.")
    return db


def get_db(db_label):
    """Get a db instance base on it's name in the config or env."""
    defaults = get_defaults()
    db_name = defaults[db_label]
    m = re.match('(\w+)://.*?/([\w.]+)', db_name)
    if m is None:
        logger.error("Poorly formed db name: %s" % db_name)
        return
    sqltype = m.groups()[0]
    return DatabaseManager(db_name, sqltype=sqltype, label=db_label)


@_clockit
def get_statements_without_agents(db, prefix, *other_stmt_clauses, **kwargs):
    """Get a generator for db orm statement objects which do not have agents."""
    num_per_yield = kwargs.pop('num_per_yield', 100)
    verbose = kwargs.pop('verbose', False)

    # Get the objects for either raw or pa statements.
    stmt_tbl_obj = db.tables[prefix + '_statements']
    agent_tbl_obj = db.tables[prefix + '_agents']

    # Build a dict mapping stmt UUIDs to statement IDs
    logger.info("Getting %s that lack %s in the database."
                % (stmt_tbl_obj.__tablename__, agent_tbl_obj.__tablename__))
    if prefix == 'pa':
        stmts_w_agents_q = db.filter_query(
            stmt_tbl_obj,
            stmt_tbl_obj.mk_hash == agent_tbl_obj.stmt_mk_hash
        )
    elif prefix == 'raw':
        stmts_w_agents_q = db.filter_query(
            stmt_tbl_obj,
            stmt_tbl_obj.id == agent_tbl_obj.stmt_id
        )
    else:
        raise IndraDatabaseError("Unrecognized prefix: %s." % prefix)
    stmts_wo_agents_q = (db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
                         .except_(stmts_w_agents_q))

    # Start printing some data
    if verbose:
        num_stmts = stmts_wo_agents_q.count()
        print("Adding agents for %d statements." % num_stmts)
    else:
        num_stmts = None

    # Get the iterator
    return stmts_wo_agents_q.yield_per(num_per_yield), num_stmts


@_clockit
def insert_agents(db, prefix, stmts_wo_agents=None, **kwargs):
    """Insert agents for statements that don't have any agents.

    Note: This method currently works for both Statements and PAStatements and
    their corresponding agents (Agents and PAAgents). However, if you already
    have preassembled INDRA Statement objects that you know don't have agents in
    the database, you can use `insert_pa_agents_directly` to insert the agents
    much faster.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are adding agents.
    prefix : str
        Select which stage of statements for which you wish to insert agents.
        The choices are 'pa' for preassembled statements or 'raw' for raw
        statements.
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
    if len(kwargs):
        raise IndraDatabaseError("Unrecognized keyword argument(s): %s."
                                 % kwargs)

    agent_tbl_obj = db.tables[prefix + '_agents']

    if stmts_wo_agents is None:
        stmts_wo_agents, num_stmts = \
            get_statements_without_agents(db, prefix, verbose=verbose)
    else:
        num_stmts = None

    if verbose:
        if num_stmts is None:
            try:
                num_stmts = len(stmts_wo_agents)
            except TypeError:
                logger.info("Could not get length from type: %s. Turning off "
                            "verbose messaging." % type(stmts_wo_agents))
                verbose = False

    # Construct the agent records
    logger.info("Building agent data for insert...")
    if verbose:
        print("Loading:", end='', flush=True)
    agent_data = []
    for i, db_stmt in enumerate(stmts_wo_agents):
        # Convert the database statement entry object into an indra statement.
        stmt = stmts_from_json([json.loads(db_stmt.json.decode())])[0]

        if prefix == 'pa':
            stmt_id = db_stmt.mk_hash
        else:  # prefix == 'raw'
            stmt_id = db_stmt.id

        agent_data.extend(_get_agent_tuples(stmt, stmt_id))

        # Optionally print another tick on the progress bar.
        if verbose and num_stmts > 25 and i % (num_stmts//25) == 0:
            print('|', end='', flush=True)

    if verbose and num_stmts > 25:
        print()

    if prefix == 'pa':
        cols = ('stmt_mk_hash', 'db_name', 'db_id', 'role')
    else:  # prefix == 'raw'
        cols = ('stmt_id', 'db_name', 'db_id', 'role')
    db.copy(agent_tbl_obj.__tablename__, agent_data, cols)
    return


@_clockit
def insert_pa_agents_directly(db, stmts, verbose=False):
    """Insert agents for preasembled statements.

    Unlike raw statements, preassembled statements are indexed by a hash,
    allowing for bulk import without a lookup beforehand, and allowing for a
    much simpler API.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are adding agents.
    stmts : list[:py:class:`Statement`]
        A list of statements for which statements should be inserted.
    verbose : bool
        If True, print extra information and a status bar while compiling
        agents for insert from statements. Default False.
    """
    if verbose:
        num_stmts = len(stmts)

    # Construct the agent records
    logger.info("Building agent data for insert...")
    if verbose:
        print("Loading:", end='', flush=True)
    agent_data = []
    for i, stmt in enumerate(stmts):
        agent_data.extend(_get_agent_tuples(stmt, stmt.get_hash(shallow=True)))

        # Optionally print another tick on the progress bar.
        if verbose and num_stmts > 25 and i % (num_stmts//25) == 0:
            print('|', end='', flush=True)

    if verbose and num_stmts > 25:
        print()

    cols = ('stmt_mk_hash', 'db_name', 'db_id', 'role')
    db.copy('pa_agents', agent_data, cols)
    return


def _get_agent_tuples(stmt, stmt_id):
    """Create the tuples for copying agents into the database."""
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
    agent_data = []
    for role, ag in agents:
        # If no agent, or no db_refs for the agent, skip the insert
        # that follows.
        if ag is None or ag.db_refs is None:
            continue
        for ns, ag_id in ag.db_refs.items():
            if isinstance(ag_id, list):
                for sub_id in ag_id:
                    agent_data.append((stmt_id, ns, sub_id, role))
            else:
                agent_data.append((stmt_id, ns, ag_id, role))
    return agent_data


def insert_db_stmts(db, stmts, db_ref_id, verbose=False):
    """Insert statement, their database, and any affiliated agents.

    Note that this method is for uploading statements that came from a
    database to our databse, not for inserting any statements to the database.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are loading statements.
    stmts : list [:py:class:`indra.statements.Statement`]
        A list of un-assembled indra statements to be uploaded to the datbase.
    db_ref_id : int
        The id to the db_ref entry corresponding to these statements.
    verbose : bool
        If True, print extra information and a status bar while compiling
        statements for insert. Default False.
    """
    # Preparing the statements for copying
    stmt_data = []
    cols = ('uuid', 'mk_hash', 'db_info_id', 'type', 'json', 'indra_version')
    if verbose:
        print("Loading:", end='', flush=True)
    for i, stmt in enumerate(stmts):
        # Only one evidence is allowed for each statement.
        for ev in stmt.evidence:
            new_stmt = stmt.make_generic_copy()
            new_stmt.evidence.append(ev)
            stmt_rec = (
                new_stmt.uuid,
                new_stmt.get_hash(),
                db_ref_id,
                new_stmt.__class__.__name__,
                json.dumps(new_stmt.to_json()).encode('utf8'),
                get_version()
            )
            stmt_data.append(stmt_rec)
        if verbose and i % (len(stmts)//25) == 0:
            print('|', end='', flush=True)
    if verbose:
        print(" Done loading %d statements." % len(stmts))
    db.copy('raw_statements', stmt_data, cols)
    stmts_to_add_agents, num_stmts = \
        get_statements_without_agents(db, 'raw',
                                      db.RawStatements.db_info_id == db_ref_id)
    insert_agents(db, 'raw', stmts_to_add_agents)
    return


def insert_pa_stmts(db, stmts, verbose=False, do_copy=True,
                    direct_agent_load=True):
    """Insert pre-assembled statements, and any affiliated agents.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are loading pre-assembled
        statements.
    stmts : iterable [:py:class:`indra.statements.Statement`]
        A list of pre-assembled indra statements to be uploaded to the datbase.
    verbose : bool
        If True, print extra information and a status bar while compiling
        statements for insert. Default False.
    direct_agent_load : bool
        If True (default), use the Statement get_hash method to get the id's of
        the Statements for insert, instead of looking up the ids of Statements
        from the database.
    """
    logger.info("Beginning to insert pre-assembled statements.")
    stmt_data = []
    indra_version = get_version()
    cols = ('uuid', 'matches_key', 'mk_hash', 'type', 'json', 'indra_version')
    if verbose:
        print("Loading:", end='', flush=True)
    for i, stmt in enumerate(stmts):
        stmt_rec = (
            stmt.uuid,
            stmt.matches_key(),
            stmt.get_hash(shallow=True),
            stmt.__class__.__name__,
            json.dumps(stmt.to_json()).encode('utf8'),
            indra_version
        )
        stmt_data.append(stmt_rec)
        if verbose and i % (len(stmts)//25) == 0:
            print('|', end='', flush=True)
    if verbose:
        print(" Done loading %d statements." % len(stmts))
    if do_copy:
        db.copy('pa_statements', stmt_data, cols)
    else:
        db.insert_many('pa_statements', stmt_data, cols=cols)
    if insert_pa_agents_directly:
        insert_pa_agents_directly(db, stmts, verbose=verbose)
    else:
        insert_agents(db, 'pa', verbose=verbose)
    return


#==============================================================================
# Below are some functions that are useful for getting raw statements from the
# database at various levels of abstraction.
#==============================================================================

def _get_statement_object(db_stmt):
    """Get an INDRA Statement object from a db_stmt."""
    return Statement._from_json(json.loads(db_stmt.json.decode('utf-8')))


def _set_evidence_text_ref(stmt, tr):
    # This is a separate function because it is likely to change, and this is a
    # critical process that is executed in multiple places.
    for ev in stmt.evidence:
        ev.pmid = tr.pmid


@_clockit
def _fix_evidence_refs(db, rid_stmt_trios):
    """Get proper id data for a raw statement from the database.

    Alterations are made to the Statement objects "in-place", so this function
    itself returns None.
    """
    rid_set = {rid for rid, _, _ in rid_stmt_trios if rid is not None}
    logger.info("Getting text refs for %d readings." % len(rid_set))
    if rid_set:
        rid_tr_pairs = db.select_all([db.Reading.id, db.TextRef],
                                     db.Reading.id.in_(rid_set),
                                     *db.join(db.TextRef, db.Reading))
        rid_tr_dict = {rid: tr for rid, tr in rid_tr_pairs}
        for rid, sid, stmt in rid_stmt_trios:
            if rid is None:
                # This means this statement came from a database, not reading.
                continue
            assert len(stmt.evidence) == 1, \
                "Only raw statements can have their refs fixed."
            _set_evidence_text_ref(stmt, rid_tr_dict[rid])
    return


@_clockit
def get_raw_stmts_frm_db_list(db, db_stmt_objs, fix_refs=True, with_sids=True):
    """Convert table objects of raw statements into INDRA Statement objects."""
    rid_stmt_sid_trios = [(db_stmt.reading_id, db_stmt.id,
                          _get_statement_object(db_stmt))
                          for db_stmt in db_stmt_objs]
    if fix_refs:
        _fix_evidence_refs(db, rid_stmt_sid_trios)
    # Note: it is important that order is maintained here (hence not a set or
    # dict).
    if with_sids:
        return [(sid, stmt) for _, sid, stmt in rid_stmt_sid_trios]
    else:
        return [stmt for _, _, stmt in rid_stmt_sid_trios]


def _get_reading_statement_dict(db, clauses=None, get_full_stmts=True):
    """Get a nested dict of statements, keyed by ref, content, and reading."""
    # Construct the query for metadata from the database.
    q = (db.session.query(db.TextRef, db.TextContent.id,
                          db.TextContent.source, db.Reading.id,
                          db.Reading.reader_version, db.RawStatements.id,
                          db.RawStatements.json)
         .filter(*db.join(db.RawStatements, db.TextRef)))
    if clauses:
        q = q.filter(*clauses)

    # Prime some counters.
    num_duplicate_evidence = 0
    num_unique_evidence = 0

    # Populate a dict with all the data.
    stmt_nd = NestedDict()
    for tr, tcid, src, rid, rv, sid, sjson in q.yield_per(1000):
        # Back out the reader name.
        for reader, rv_list in reader_versions.items():
            if rv in rv_list:
                break
        else:
            raise Exception("rv %s not recognized." % rv)

        # Get the json for comparison and/or storage
        stmt_json = json.loads(sjson.decode('utf8'))
        stmt = Statement._from_json(stmt_json)
        _set_evidence_text_ref(stmt, tr)

        # Hash the compbined stmt and evidence matches key.
        stmt_hash = stmt.get_hash()

        # For convenience get the endpoint statement dict
        s_dict = stmt_nd[tr.id][src][tcid][reader][rv][rid]

        # Initialize the value to a set, and count duplicates
        if stmt_hash not in s_dict.keys():
            s_dict[stmt_hash] = set()
            num_unique_evidence += 1
        else:
            num_duplicate_evidence += 1

        # Either store the statement, or the statement id.
        if get_full_stmts:
            s_dict[stmt_hash].add((sid, stmt))
        else:
            s_dict[stmt_hash].add((sid, None))

    # Report on the results.
    print("Found %d relevant text refs with statements." % len(stmt_nd))
    print("number of statement exact duplicates: %d" % num_duplicate_evidence)
    print("number of unique statements: %d" % num_unique_evidence)
    return stmt_nd


# Specify versions of readers, and preference. Later in the list is better.
reader_versions = {
    'sparser': ['sept14-linux\n', 'sept14-linux', 'June2018-linux'],
    'reach': ['61059a-biores-e9ee36', '1.3.3-61059a-biores-']
    }

# Specify sources of fulltext content, and order priorities.
text_content_sources = ['pubmed', 'elsevier', 'manuscripts', 'pmc_oa']


def _get_filtered_rdg_statements(stmt_nd, get_full_stmts, linked_sids=None,
                                 ignore_duplicates=False):
    """Get the set of statements/ids from readings minus exact duplicates."""
    logger.info("Filtering the statements from reading.")
    if linked_sids is None:
        linked_sids = set()
    def better_func(element):
        return text_content_sources.index(element)

    # Now we filter and get the set of statements/statement ids.
    stmt_tpls = set()
    duplicate_sids = set()  # Statements that are exact duplicates.
    bettered_duplicate_sids = set()  # Statements with "better" alternatives
    for trid, src_dict in stmt_nd.items():
        some_bettered_duplicate_tpls = set()
        # Filter out unneeded fulltext.
        while len(src_dict) > 1:
            try:
                worst_src = min(src_dict, key=better_func)
                some_bettered_duplicate_tpls |= src_dict[worst_src].get_leaves()
                del src_dict[worst_src]
            except:
                print(src_dict)
                raise

        # Filter out the older reader versions
        for reader, rv_list in reader_versions.items():
            for rv_dict in src_dict.gets(reader):
                best_rv = max(rv_dict, key=lambda x: rv_list.index(x))

                # Record the rest of the statement uuids.
                for rv, r_dict in rv_dict.items():
                    if rv != best_rv:
                        some_bettered_duplicate_tpls |= r_dict.get_leaves()

                # Take any one of the duplicates. Statements/Statement ids are
                # already grouped into sets of duplicates keyed by the
                # Statement and Evidence matches key hashes. We only want one
                # of each.
                stmt_set_itr = (stmt_set for r_dict in rv_dict[best_rv].values()
                                for stmt_set in r_dict.values())
                if ignore_duplicates:
                    some_stmt_tpls = {stmt_tpl for stmt_set in stmt_set_itr
                                      for stmt_tpl in stmt_set}
                else:
                    some_stmt_tpls, some_duplicate_tpls = \
                        _detect_exact_duplicates(stmt_set_itr, linked_sids)

                    # Get the sids for the statements.
                    duplicate_sids |= {sid for sid, _ in some_duplicate_tpls}

                stmt_tpls |= some_stmt_tpls

        # Add the bettered duplicates found in this round.
        bettered_duplicate_sids |= \
            {sid for sid, _ in some_bettered_duplicate_tpls}

    if get_full_stmts:
        stmts = {stmt for _, stmt in stmt_tpls if stmt is not None}
        assert len(stmts) == len(stmt_tpls),\
            ("Some statements were None! The interaction between "
             "_get_reading_statement_dict and _filter_rdg_statements was "
             "probably mishandled.")
    else:
        stmts = {sid for sid, _ in stmt_tpls}

    return stmts, duplicate_sids, bettered_duplicate_sids


def _detect_exact_duplicates(stmt_set_itr, linked_sids):
    # Pick one among any exact duplicates. Unlike with bettered
    # duplicates, these choices are arbitrary, and such duplicates
    # can be deleted.
    stmt_tpls = set()
    some_duplicate_tpls = set()
    for stmt_tpl_set in stmt_set_itr:
        if not stmt_tpl_set:
            continue
        elif len(stmt_tpl_set) == 1:
            # There isn't really a choice here.
            stmt_tpls |= stmt_tpl_set
        else:
            prefed_tpls = {tpl for tpl in stmt_tpl_set
                           if tpl[0] in linked_sids}
            if not prefed_tpls:
                # Pick the first one to pop, record the rest as
                # duplicates.
                stmt_tpls.add(stmt_tpl_set.pop())
                some_duplicate_tpls |= stmt_tpl_set
            elif len(prefed_tpls) == 1:
                # There is now no choice: just take the preferred
                # statement.
                stmt_tpls |= prefed_tpls
                some_duplicate_tpls |= (stmt_tpl_set - prefed_tpls)
            else:
                # This shouldn't happen, so an early run of this
                # function must have failed somehow, or else there
                # was some kind of misuse. Flag it, pick just one of
                # the preferred statements, and delete any deletable
                # statements.
                assert False, \
                    ("Duplicate deduplicated statements found: %s"
                     % str(prefed_tpls))
    return stmt_tpls, some_duplicate_tpls


def _choose_unique(not_duplicates, get_full_stmts, stmt_tpl_grp):
    """Choose one of the statements from a redundant set."""
    assert stmt_tpl_grp, "This cannot be empty."
    if len(stmt_tpl_grp) == 1:
        s_tpl = stmt_tpl_grp[0]
        duplicate_ids = set()
    else:
        stmt_tpl_set = set(stmt_tpl_grp)
        preferred_tpls = {tpl for tpl in stmt_tpl_set
                          if tpl[1] in not_duplicates}
        if not preferred_tpls:
            s_tpl = stmt_tpl_set.pop()
        elif len(preferred_tpls) == 1:
            s_tpl = preferred_tpls.pop()
        else:  # len(preferred_stmts) > 1
            assert False, \
                ("Duplicate deduplicated statements found: %s"
                 % str(preferred_tpls))
        duplicate_ids = {tpl[1] for tpl in stmt_tpl_set
                         if tpl[1] not in not_duplicates}

    if get_full_stmts:
        stmt_json = json.loads(s_tpl[2].decode('utf-8'))
        ret_stmt = Statement._from_json(stmt_json)
    else:
        ret_stmt = s_tpl[1]
    return ret_stmt, duplicate_ids


def _get_filtered_db_statements(db, get_full_stmts=False, clauses=None,
                                not_duplicates=None, num_procs=1):
    """Get the set of statements/ids from databases minus exact duplicates."""
    if not_duplicates is None:
        not_duplicates = set()

    # Only get the json if it's going to be used. Note that if the use of the
    # get_full_stmts parameter is inconsistent in _choose_unique, this will
    # cause some problems.
    if get_full_stmts:
        tbl_list = [db.RawStatements.mk_hash, db.RawStatements.id,
                    db.RawStatements.json]
    else:
        tbl_list = [db.RawStatements.mk_hash, db.RawStatements.id]

    db_s_q = db.filter_query(tbl_list, db.RawStatements.db_info_id.isnot(None))

    # Add any other criterion specified at higher levels.
    if clauses:
        db_s_q = db_s_q.filter(*clauses)

    # Produce a generator of statement groups.
    db_stmt_data = db_s_q.order_by(db.RawStatements.mk_hash).yield_per(10000)
    choose_unique_stmt = partial(_choose_unique, not_duplicates, get_full_stmts)
    stmt_groups = (list(grp) for _, grp
                   in groupby(db_stmt_data, key=lambda x: x[0]))

    # Actually do the comparison.
    if num_procs is 1:
        stmts = set()
        duplicate_ids = set()
        for stmt_list in stmt_groups:
            stmt, some_duplicates = choose_unique_stmt(stmt_list)
            stmts.add(stmt)
            duplicate_ids |= some_duplicates
    else:
        pool = Pool(num_procs)
        print("Filtering db statements in %d processess." % num_procs)
        res = pool.map(choose_unique_stmt, stmt_groups)
        pool.close()
        pool.join()
        stmt_list, duplicate_sets = zip(*res)
        stmts = set(stmt_list)
        duplicate_ids = {uuid for dup_set in duplicate_sets for uuid in dup_set}

    return stmts, duplicate_ids


@_clockit
def distill_stmts(db, get_full_stmts=False, clauses=None, num_procs=1,
                  delete_duplicates=True, weed_evidence=True, batch_size=1000):
    """Get a corpus of statements from clauses and filters duplicate evidence.

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        A database manager instance to access the database.
    get_full_stmts : bool
        By default (False), only Statement ids (the primary index of Statements
        on the database) are returned. However, if set to True, serialized
        INDRA Statements will be returned. Note that this will in general be
        VERY large in memory, and therefore should be used with caution.
    clauses : None or list of sqlalchemy clauses
        By default None. Specify sqlalchemy clauses to reduce the scope of
        statements, e.g. `clauses=[db.Statements.type == 'Phosphorylation']` or
        `clauses=[db.Statements.uuid.in_([<uuids>])]`.
    num_procs : int
        Select the number of process that can be used.
    delete_duplicates : bool
        Choose whether you want to delete the statements that are found to be
        duplicates.
    weed_evidence : bool
        If True, evidence links that exist for raw statements that now have
        better alternatives will be removed. If False, such links will remain,
        which may cause problems in incremental pre-assembly.

    Returns
    -------
    stmt_ret : set
        A set of either statement ids or serialized statements, depending on
        `get_full_stmts`.
    """
    if delete_duplicates:
        logger.info("Looking for ids from existing links...")
        linked_sids = {sid for sid,
                        in db.select_all(db.RawUniqueLinks.raw_stmt_id)}
    else:
        linked_sids = set()

    # Get de-duplicated Statements, and duplicate uuids, as well as uuid of
    # Statements that have been improved upon...
    logger.info("Sorting reading statements...")
    stmt_nd = _get_reading_statement_dict(db, clauses, get_full_stmts)

    stmts, duplicate_sids, bettered_duplicate_sids = \
        _get_filtered_rdg_statements(stmt_nd, get_full_stmts, linked_sids)
    logger.info("After filtering reading: %d unique statements, %d exact "
                "duplicates, %d with results from better resources available."
                % (len(stmts), len(duplicate_sids),
                   len(bettered_duplicate_sids)))
    assert not linked_sids & duplicate_sids, linked_sids & duplicate_sids
    del stmt_nd  # This takes up a lot of memory, and is done being used.

    db_stmts, db_duplicates = \
        _get_filtered_db_statements(db, get_full_stmts, clauses, linked_sids,
                                    num_procs)
    stmts |= db_stmts
    duplicate_sids |= db_duplicates
    logger.info("After filtering database statements: %d unique, %d duplicates."
                % (len(stmts), len(duplicate_sids)))
    assert not linked_sids & duplicate_sids, linked_sids & duplicate_sids

    # Remove support links for statements that have better versions available.
    bad_link_sids = bettered_duplicate_sids & linked_sids
    if len(bad_link_sids) and weed_evidence:
        logger.info("Removing bettered evidence links...")
        rm_links = db.select_all(
            db.RawUniqueLinks,
            db.RawUniqueLinks.raw_stmt_id.in_(bad_link_sids)
            )
        db.delete_all(rm_links)

    # Delete exact duplicates
    if len(duplicate_sids) and delete_duplicates:
        logger.info("Deleting duplicates...")
        for dup_id_batch in batch_iter(duplicate_sids, batch_size, set):
            bad_stmts = db.select_all(db.RawStatements,
                                      db.RawStatements.id.in_(dup_id_batch))
            bad_sid_set = {s.id for s in bad_stmts}
            bad_agents = db.select_all(db.RawAgents,
                                       db.RawAgents.stmt_id.in_(bad_sid_set))
            logger.info("Deleting %d agents associated with redundant raw "
                        "statements." % len(bad_agents))
            db.delete_all(bad_agents)
            logger.info("Deleting %d redundant raw statements." % len(bad_stmts))
            db.delete_all(bad_stmts)

    return stmts


def unpack(bts, decode=True):
    ret = zlib.decompress(bts, zlib.MAX_WBITS+16)
    if decode:
        ret = ret.decode('utf-8')
    return ret
