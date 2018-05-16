from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['get_defaults', 'get_primary_db', 'insert_agents', 'insert_pa_stmts',
           'insert_db_stmts', 'make_stmts_from_db_list']

import re
import json
import logging
import zlib
from indra.util.get_version import get_version
from indra.statements import Complex, SelfModification, ActiveForm,\
    stmts_from_json, Conversion, Translocation, Statement
from .database_manager import DatabaseManager, IndraDatabaseError, texttypes
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
    for k in key_list:
        test_name = test_defaults[k]
        m = re.match('(\w+)://.*?/([\w.]+)', test_name)
        sqltype = m.groups()[0]
        try:
            db = DatabaseManager(test_name, sqltype=sqltype)
            db.grab_session()
        except Exception as e:
            logger.error("%s didn't work" % test_name)
            logger.exception(e)
            continue  # Clearly this test database won't work.
        logger.info("Using test database %s." % k)
        break
    else:
        logger.error("Could not load a test database!")
        return None
    return db


def insert_agents(db, stmt_tbl_obj, agent_tbl_obj, *other_stmt_clauses,
                  **kwargs):
    """Insert agents for statements that don't have any agents.

    Note: This method currently works for both Statements and PAStatements and
    their corresponding agents (Agents and PAAgents).

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are adding agents.
    stmt_tbl_obj : :py:class:`sqlalchemy.Base` table object
        For example, `db.Statements`. The object corresponding to the
        statements column you creating agents for.
    agent_tbl_obj : :py:class:`sqlalchemy.Base` table object
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
    override_default_query = kwargs.pop('override_default_query', False)
    if len(kwargs):
        raise IndraDatabaseError("Unrecognized keyword argument(s): %s."
                                 % kwargs)
    # Build a dict mapping stmt UUIDs to statement IDs
    logger.info("Getting %s that lack %s in the database."
                % (stmt_tbl_obj.__tablename__, agent_tbl_obj.__tablename__))
    if not override_default_query:
        stmts_w_agents_q = db.filter_query(
            stmt_tbl_obj,
            stmt_tbl_obj.id == agent_tbl_obj.stmt_id
            )
        stmts_wo_agents_q = (db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
                             .except_(stmts_w_agents_q))
    else:
        stmts_wo_agents_q = db.filter_query(stmt_tbl_obj, *other_stmt_clauses)
    logger.debug("Getting stmts with query:\n%s" % str(stmts_wo_agents_q))
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
                if isinstance(ag_id, list):
                    for sub_id in ag_id:
                        agent_data.append((db_stmt.id, ns, sub_id, role))
                else:
                    agent_data.append((db_stmt.id, ns, ag_id, role))

        # Optionally print another tick on the progress bar.
        if verbose and num_stmts > 25 and i % (num_stmts//25) == 0:
            print('|', end='', flush=True)

    if verbose and num_stmts > 25:
        print()

    cols = ('stmt_id', 'db_name', 'db_id', 'role')
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

    Parameters
    ----------
    db : :py:class:`DatabaseManager`
        The manager for the database into which you are loading pre-assembled
        statements.
    stmts : list [:py:class:`indra.statements.Statement`]
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


def make_stmts_from_db_list(db_stmt_objs):
    stmt_json_list = []
    for st_obj in db_stmt_objs:
        stmt_json_list.append(json.loads(st_obj.json.decode('utf8')))
    return stmts_from_json(stmt_json_list)


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


def distill_stmts_from_reading(db, get_full_stmts=False, clauses=None):
    """Get a corpus of statements from clauses and filters duplicate evidence.

    Note that this will only get statements from reading.

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

    Returns
    -------
    stmt_dn : NestedDict
        A deeply nested recursive dictionary, carrying the metadata for the
        Statements.
    stmt_ret : set
        A set of either statement ids or serialized statements, depending on
        `get_full_stmts`.
    """
    # Construct the query for metadata from the database.
    q = (db.session.query(db.TextContent.text_ref_id, db.TextContent.id,
                          db.TextContent.source, db.Readings.id,
                          db.Readings.reader_version, db.Statements.id,
                          db.Statements.json)
         .filter(db.TextContent.id == db.Readings.text_content_id,
                 db.Readings.id == db.Statements.reader_ref))
    if clauses:
        q.filter(*clauses)

    # Specify sources of fulltext content, and order priorities.
    full_text_content = ['manuscripts', 'pmc_oa', 'elsevier']

    # Specify versions of readers, and preference.
    sparser_versions = ['sept14-linux\n', 'sept14-linux']
    reach_versions = ['61059a-biores-e9ee36', '1.3.3-61059a-biores-']

    # Prime some counters.
    num_duplicate_evidence = 0
    num_unique_evidence = 0

    # Populate a dict with all the data.
    stmt_nd = NestedDict()
    for trid, tcid, src, rid, rv, sid, sjson in q.yield_per(1000):
        # Back out the reader name.
        if rv in sparser_versions:
            reader = 'sparser'
        elif rv in reach_versions:
            reader = 'reach'
        else:
            raise Exception("rv %s not recognized." % rv)

        # Get the json for comparison and/or storage
        stmt_json = json.loads(sjson.decode('utf8'))
        stmt = Statement._from_json(stmt_json)

        # Hash the compbined stmt and evidence matches key.
        m_key = stmt.matches_key() + stmt.evidence[0].matches_key()
        stmt_hash = hash(m_key)

        # For convenience get the endpoint statement dict
        s_dict = stmt_nd[trid][src][tcid][reader][rv][rid]

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

    # Now we filter and get the set of statements/statement ids.
    stmts = set()
    for trid, src_dict in stmt_nd.items():
        # Filter out unneeded fulltext.
        while sum([k != 'pubmed' for k in src_dict.keys()]) > 1:
            worst_src = min(src_dict,
                            key=lambda x: full_text_content.index(x[0]))
            del src_dict[worst_src]

        # Filter out the older reader versions
        for reader, rv_list in [('reach', reach_versions),
                                ('sparser', sparser_versions)]:
            for rv_dict in src_dict.gets(reader):
                best_rv = max(rv_dict, key=lambda x: rv_list.index(x))

                # Take any one of the duplicates. Statements/Statement ids are
                # already grouped into sets of duplicates keyed by the
                # Statement and Evidence matches key hashes. We only want one
                # of each.
                stmts |= {(ev_hash, list(ev_set)[0])
                          for ev_hash, ev_set in rv_dict[best_rv].items()}

    return stmt_nd, stmts


def unpack(bts, decode=True):
    ret = zlib.decompress(bts, zlib.MAX_WBITS+16)
    if decode:
        ret = ret.decode('utf-8')
    return ret
