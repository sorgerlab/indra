from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str


__all__ = ['get_defaults', 'get_primary_db', 'get_db', 'insert_agents',
           'insert_pa_stmts', 'insert_db_stmts', 'make_raw_stmts_from_db_list',
           'distill_stmts']

import re
import json
import zlib
import logging
from itertools import groupby
from functools import partial
from multiprocessing.pool import Pool

from indra.util import batch_iter
from indra.util.get_version import get_version
from indra.statements import Complex, SelfModification, ActiveForm,\
    stmts_from_json, Conversion, Translocation, Statement

from .database_manager import DatabaseManager, IndraDatabaseError
from .config import get_databases as get_defaults

logger = logging.getLogger('db_util')


__PRIMARY_DB = None


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


def insert_agents(db, prefix, *other_stmt_clauses, **kwargs):
    """Insert agents for statements that don't have any agents.

    Note: This method currently works for both Statements and PAStatements and
    their corresponding agents (Agents and PAAgents).

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
    stmt_tbl_obj = db.tables[prefix + '_statements']
    agent_tbl_obj = db.tables[prefix + '_agents']
    num_per_yield = kwargs.pop('num_per_yield', 100)
    override_default_query = kwargs.pop('override_default_query', False)
    if len(kwargs):
        raise IndraDatabaseError("Unrecognized keyword argument(s): %s."
                                 % kwargs)
    # Build a dict mapping stmt UUIDs to statement IDs
    logger.info("Getting %s that lack %s in the database."
                % (stmt_tbl_obj.__tablename__, agent_tbl_obj.__tablename__))
    if not override_default_query:
        if prefix == 'pa':
            stmts_w_agents_q = db.filter_query(
                stmt_tbl_obj,
                stmt_tbl_obj.mk_hash == agent_tbl_obj.stmt_mk_hash
            )
        elif prefix == 'raw':
            stmts_w_agents_q = db.filter_query(
                stmt_tbl_obj,
                stmt_tbl_obj.uuid == agent_tbl_obj.stmt_uuid
                )
        stmts_wo_agents_q = (db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
                             .except_(stmts_w_agents_q))
    else:
        stmts_wo_agents_q = db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
    #logger.debug("Getting stmts with query:\n%s" % str(stmts_wo_agents_q))
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
        stmt = stmts_from_json([json.loads(db_stmt.json.decode())])[0]

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
                if prefix == 'pa':
                    stmt_id = db_stmt.mk_hash
                elif prefix == 'raw':
                    stmt_id = db_stmt.uuid
                if isinstance(ag_id, list):
                    for sub_id in ag_id:
                        agent_data.append((stmt_id, ns, sub_id, role))
                else:
                    agent_data.append((stmt_id, ns, ag_id, role))

        # Optionally print another tick on the progress bar.
        if verbose and num_stmts > 25 and i % (num_stmts//25) == 0:
            print('|', end='', flush=True)

    if verbose and num_stmts > 25:
        print()

    if prefix == 'pa':
        cols = ('stmt_mk_hash', 'db_name', 'db_id', 'role')
    elif prefix == 'raw':
        cols = ('stmt_uuid', 'db_name', 'db_id', 'role')
    db.copy(agent_tbl_obj.__tablename__, agent_data, cols)
    return


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
    insert_agents(db, 'raw', db.RawStatements.db_info_id == db_ref_id)
    return


def insert_pa_stmts(db, stmts, verbose=False, do_copy=True):
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
    """
    logger.info("Beginning to insert pre-assembled statements.")
    stmt_data = []
    indra_version = get_version()
    cols = ('uuid', 'mk_hash', 'type', 'json', 'indra_version')
    if verbose:
        print("Loading:", end='', flush=True)
    for i, stmt in enumerate(stmts):
        stmt_rec = (
            stmt.uuid,
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


def _fix_evidence_refs(db, rid_stmt_pairs):
    """Get proper id data for a raw statement from the database.

    Alterations are made to the Statement objects "in-place", so this function
    itself returns None.
    """
    rid_set = {rid for rid, _ in rid_stmt_pairs if rid is not None}
    logger.info("Getting text refs for %d readings." % len(rid_set))
    if rid_set:
        rid_tr_pairs = db.select_all([db.Reading.id, db.TextRef],
                                     db.Reading.id.in_(rid_set),
                                     *db.join(db.TextRef, db.Reading))
        rid_tr_dict = {rid: tr for rid, tr in rid_tr_pairs}
        for rid, stmt in rid_stmt_pairs:
            if rid is None:
                # This means this statement came from a database, not reading.
                continue
            assert len(stmt.evidence) == 1, \
                "Only raw statements can have their refs fixed."
            _set_evidence_text_ref(stmt, rid_tr_dict[rid])
    return


def make_raw_stmts_from_db_list(db, db_stmt_objs):
    """Convert table objects of raw statements into INDRA Statement objects."""
    rid_stmt_pairs = [(db_stmt.reading_id, _get_statement_object(db_stmt))
                      for db_stmt in db_stmt_objs]
    _fix_evidence_refs(db, rid_stmt_pairs)
    return [stmt for _, stmt in rid_stmt_pairs]


class NestedDict(dict):
    """A dict-like object that recursively populates elements of a dict."""

    def __getitem__(self, key):
        if key not in self.keys():
            val = self.__class__()
            self.__setitem__(key, val)
        else:
            val = dict.__getitem__(self, key)
        return val

    def __repr__(self):
        sub_str = dict.__repr__(self)[1:-1]
        if not sub_str:
            return self.__class__.__name__ + '()'
        # This does not completely generalize, but it works for most cases.
        for old, new in [('), ', '),\n'), ('\n', '\n  ')]:
            sub_str = sub_str.replace(old, new)
        return'%s(\n  %s\n)' % (self.__class__.__name__, sub_str)

    def __str__(self):
        return self.__repr__()

    def export_dict(self):
        "Convert this into an ordinary dict (of dicts)."
        return {k: v.export_dict() if isinstance(v, self.__class__) else v
                for k, v in self.items()}

    def get(self, key):
        "Find the first value within the tree which has the key."
        if key in self.keys():
            return self[key]
        else:
            res = None
            for v in self.values():
                # This could get weird if the actual expected returned value
                # is None, especially in teh case of overlap. Any ambiguity
                # would be resolved by get_path(s).
                if hasattr(v, 'get'):
                    res = v.get(key)
                if res is not None:
                    break
            return res

    def get_path(self, key):
        "Like `get`, but also return the path taken to the value."
        if key in self.keys():
            return (key,), self[key]
        else:
            key_path, res = (None, None)
            for sub_key, v in self.items():
                if isinstance(v, self.__class__):
                    key_path, res = v.get_path(key)
                elif hasattr(v, 'get'):
                    res = v.get(key)
                    key_path = (key,) if res is not None else None
                if res is not None and key_path is not None:
                    key_path = (sub_key,) + key_path
                    break
            return key_path, res

    def gets(self, key):
        "Like `get`, but return all matches, not just the first."
        result_list = []
        if key in self.keys():
            result_list.append(self[key])
        for v in self.values():
            if isinstance(v, self.__class__):
                sub_res_list = v.gets(key)
                for res in sub_res_list:
                    result_list.append(res)
            elif isinstance(v, dict):
                if key in v.keys():
                    result_list.append(v[key])
        return result_list

    def get_paths(self, key):
        "Like `gets`, but include the paths, like `get_path` for all matches."
        result_list = []
        if key in self.keys():
            result_list.append(((key,), self[key]))
        for sub_key, v in self.items():
            if isinstance(v, self.__class__):
                sub_res_list = v.get_paths(key)
                for key_path, res in sub_res_list:
                    result_list.append(((sub_key,) + key_path, res))
            elif isinstance(v, dict):
                if key in v.keys():
                    result_list.append(((sub_key, key), v[key]))
        return result_list

    def get_leaves(self):
        """Get the deepest entries as a flat set."""
        ret_set = set()
        for val in self.values():
            if isinstance(val, self.__class__):
                ret_set |= val.get_leaves()
            elif isinstance(val, dict):
                ret_set |= set(val.values())
            elif isinstance(val, list):
                ret_set |= set(val)
            elif isinstance(val, set):
                ret_set |= val
            else:
                ret_set.add(val)
        return ret_set


def _get_reading_statement_dict(db, reader_versions, get_full_stmts=False,
                                clauses=None):
    """Get a nested dict of statements, keyed by ref, content, and reading."""
    # Construct the query for metadata from the database.
    q = (db.session.query(db.TextRef, db.TextContent.id,
                          db.TextContent.source, db.Reading.id,
                          db.Reading.reader_version, db.RawStatements.uuid,
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
            s_dict[stmt_hash].add(stmt)
        else:
            s_dict[stmt_hash].add(sid)

    # Report on the results.
    print("Found %d relevant text refs with statements." % len(stmt_nd))
    print("number of statement exact duplicates: %d" % num_duplicate_evidence)
    print("number of unique statements: %d" % num_unique_evidence)
    return stmt_nd


def _get_filtered_rdg_statements(db, get_full_stmts, clauses=None,
                                 linked_uuids=None):
    """Get the set of statements/ids from readings minus exact duplicates."""
    if linked_uuids is None:
        linked_uuids = set()
    # Specify versions of readers, and preference.
    reader_versions = {
        'sparser': ['sept14-linux\n', 'sept14-linux'],
        'reach': ['61059a-biores-e9ee36', '1.3.3-61059a-biores-']
    }

    stmt_nd = _get_reading_statement_dict(db, reader_versions, get_full_stmts,
                                          clauses)

    # Specify sources of fulltext content, and order priorities.
    full_text_content = ['manuscripts', 'pmc_oa', 'elsevier', 'pubmed']

    def better_func(element):
        return full_text_content.index(element)

    # Now we filter and get the set of statements/statement ids.
    stmts = set()
    duplicate_uuids = set()  # Statements that are exact duplicates.
    bettered_duplicate_uuids = set()  # Statements with "better" alternatives
    for trid, src_dict in stmt_nd.items():
        bettered_duplicate_stmts = set()
        # Filter out unneeded fulltext.
        while len(src_dict) > 1:
            try:
                worst_src = min(src_dict, key=better_func)
                bettered_duplicate_stmts |= src_dict[worst_src].get_leaves()
                del src_dict[worst_src]
            except:
                print(src_dict)
                raise

        # Filter out the older reader versions
        for reader, rv_list in reader_versions.items():
            for rv_dict in src_dict.gets(reader):
                best_rv = max(rv_dict, key=lambda x: rv_list.index(x))

                # Take any one of the duplicates. Statements/Statement ids are
                # already grouped into sets of duplicates keyed by the
                # Statement and Evidence matches key hashes. We only want one
                # of each.
                stmt_set_itr = (stmt_set for r_dict in rv_dict[best_rv].values()
                                for stmt_set in r_dict.values())

                # Record the rest of the statement uuids.
                for rv, r_dict in rv_dict.items():
                    if rv != best_rv:
                        bettered_duplicate_stmts |= r_dict.get_leaves()

                # Pick one among any exact duplicates. Unlike with bettered
                # duplicates, these choices are arbitrary, and such duplicates
                # can be deleted.
                duplicate_stmts = set()
                for stmt_set in stmt_set_itr:
                    if not stmt_set:
                        continue
                    elif len(stmt_set) == 1:
                        # There isn't really a choice here.
                        stmts |= stmt_set
                    else:
                        if get_full_stmts:
                            preferred_stmts = {s for s in stmt_set
                                               if s.uuid in linked_uuids}
                        else:
                            preferred_stmts = stmt_set & linked_uuids
                        if not preferred_stmts:
                            # Pick the first one to pop, record the rest as
                            # duplicates.
                            stmts.add(stmt_set.pop())
                            duplicate_stmts |= stmt_set
                        elif len(preferred_stmts) == 1:
                            # There is now no choice: just take the preferred
                            # statement.
                            stmts |= preferred_stmts
                            duplicate_stmts |= (stmt_set - preferred_stmts)
                        else:
                            # This shouldn't happen, so an early run of this
                            # function must have failed somehow, or else there
                            # was some kind of misuse. Flag it, pick just one of
                            # the preferred statements, and delete any deletable
                            # statements.
                            assert False,\
                                ("Duplicate deduplicated statements found: %s"
                                 % str(preferred_stmts))

                # Get the uuids from the statements, if we have full statements.
                if get_full_stmts:
                    duplicate_uuids |= {s.uuid for s in duplicate_stmts}
                    bettered_duplicate_uuids |= \
                        {s.uuid for s in bettered_duplicate_stmts}
                else:
                    duplicate_uuids |= duplicate_stmts
                    bettered_duplicate_uuids |= bettered_duplicate_stmts

    return stmts, duplicate_uuids, bettered_duplicate_uuids


def _choose_unique(not_duplicates, get_full_stmts, stmt_tpl_grp):
    """Choose one of the statements from a redundant set."""
    assert stmt_tpl_grp, "This cannot be empty."
    if len(stmt_tpl_grp) == 1:
        s_tpl = stmt_tpl_grp[0]
        duplicate_ids = set()
    else:
        stmt_tpl_set = set(stmt_tpl_grp)
        preferred_stmts = {s for s in stmt_tpl_set if s.uuid in not_duplicates}
        if not preferred_stmts:
            s_tpl = stmt_tpl_set.pop()
        elif len(preferred_stmts) == 1:
            s_tpl = preferred_stmts.pop()
        else:  # len(preferred_stmts) > 1
            assert False, \
                ("Duplicate deduplicated statements found: %s"
                 % str(preferred_stmts))
        duplicate_ids = {s.uuid for s in stmt_tpl_set
                         if s.uuid not in not_duplicates}

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

    db_s_q = db.filter_query([db.RawStatements.mk_hash,
                              db.RawStatements.uuid,
                              db.RawStatements.json],
                             db.RawStatements.db_info_id.isnot(None))
    if clauses:
        db_s_q = db_s_q.filter(*clauses)
    db_stmt_data = db_s_q.order_by(db.RawStatements.mk_hash).yield_per(10000)
    choose_unique_stmt = partial(_choose_unique, not_duplicates, get_full_stmts)
    stmt_groups = (list(grp) for _, grp
                   in groupby(db_stmt_data, key=lambda x: x[0]))
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
        linked_uuids = {uuid for uuid,
                        in db.select_all(db.RawUniqueLinks.raw_stmt_uuid)}
    else:
        linked_uuids = set()

    # Get de-duplicated Statements, and duplicate uuids, as well as uuid of
    # Statements that have been improved upon...
    stmts, duplicate_uuids, bettered_duplicate_uuids = \
        _get_filtered_rdg_statements(db, get_full_stmts, clauses, linked_uuids)
    print("After filtering reading: %d unique statements, %d duplicates."
          % (len(stmts), len(duplicate_uuids)))
    assert not linked_uuids & duplicate_uuids, linked_uuids & duplicate_uuids

    db_stmts, db_duplicates = \
        _get_filtered_db_statements(db, get_full_stmts, clauses, linked_uuids,
                                    num_procs)
    stmts |= db_stmts
    duplicate_uuids |= db_duplicates
    print("After filtering database statements: %d unique, %d duplicates."
          % (len(stmts), len(duplicate_uuids)))
    assert not linked_uuids & duplicate_uuids, linked_uuids & duplicate_uuids

    # Remove support links for statements that have better versions available.
    bad_link_uuids = bettered_duplicate_uuids & linked_uuids
    if len(bad_link_uuids) and weed_evidence:
        print("Removing bettered evidence links...")
        rm_links = db.select_all(
            db.RawUniqueLinks,
            db.RawUniqueLinks.raw_stmt_uuid.in_(bad_link_uuids)
            )
        db.delete_all(rm_links)

    # Delete exact duplicates
    if len(duplicate_uuids) and delete_duplicates:
        print("Deleting duplicates...")
        for dup_id_batch in batch_iter(duplicate_uuids, batch_size, set):
            bad_stmts = db.select_all(db.RawStatements,
                                      db.RawStatements.uuid.in_(dup_id_batch))
            bad_uuid_set = {s.uuid for s in bad_stmts}
            bad_agents = db.select_all(db.RawAgents,
                                       db.RawAgents.stmt_uuid.in_(bad_uuid_set))
            print("Deleting %d agents associated with redundant raw statements."
                  % len(bad_agents))
            db.delete_all(bad_agents)
            print("Deleting %d redundant raw statements." % len(bad_stmts))
            db.delete_all(bad_stmts)

    return stmts


def unpack(bts, decode=True):
    ret = zlib.decompress(bts, zlib.MAX_WBITS+16)
    if decode:
        ret = ret.decode('utf-8')
    return ret
