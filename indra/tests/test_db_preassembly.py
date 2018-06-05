from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
import pickle
import random
import logging
from datetime import timedelta, datetime


gm_logger = logging.getLogger('grounding_mapper')
gm_logger.setLevel(logging.WARNING)

sm_logger = logging.getLogger('sitemapper')
sm_logger.setLevel(logging.WARNING)

ps_logger = logging.getLogger('phosphosite')
ps_logger.setLevel(logging.WARNING)

pa_logger = logging.getLogger('preassembler')
pa_logger.setLevel(logging.WARNING)

from indra.db import util as db_util
from indra.db import client as db_client
from indra.db import preassembly_manager as pm
from indra.statements import stmts_from_json, Statement
from indra.tools import assemble_corpus as ac

from nose.plugins.attrib import attr
from .util import needs_py3

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MAX_NUM_STMTS = 11721
BATCH_SIZE = 2017
STMTS = None


def _load_tuples(fname):
    with open(os.path.join(THIS_DIR, fname), 'rb') as f:
        ret_tuples = pickle.load(f)
    return ret_tuples


def _get_stmt_tuples():
    stmt_tuples = _load_tuples('test_stmts_tuples.pkl')
    col_names = stmt_tuples.pop(0)
    if len(stmt_tuples) > MAX_NUM_STMTS:
        stmt_tuples = random.sample(stmt_tuples, MAX_NUM_STMTS)
    return stmt_tuples, col_names


def _get_background_loaded_db():
    db = db_util.get_test_db()
    db._clear(force=True)

    # Get and load the provenance for the statements.
    print("\tLoading background metadata...")
    db.copy('text_ref', _load_tuples('test_text_ref_tuples.pkl'),
            ('id', 'pmid', 'pmcid', 'doi'))
    tc_tuples = [t + (b'',)
                 for t in _load_tuples('test_text_content_tuples.pkl')]
    db.copy('text_content', tc_tuples, ('id', 'text_ref_id', 'source', 'format',
                                        'text_type', 'content'))
    r_tuples = [t + (b'',) for t in _load_tuples('test_reading_tuples.pkl')]
    db.copy('reading', r_tuples, ('id', 'reader', 'reader_version',
                                  'text_content_id', 'format', 'bytes'))
    db.copy('db_info', _load_tuples('test_db_info_tuples.pkl'),
            ('id', 'db_name'))
    return db


def _get_input_stmt_tuples():
    print("\tPrepping the raw statements...")
    stmt_tuples, col_names = _get_stmt_tuples()
    copy_col_names = ('uuid', 'mk_hash', 'type', 'indra_version', 'json',
                      'reading_id', 'db_info_id')
    copy_stmt_tuples = []
    for tpl in stmt_tuples:
        entry_dict = dict(zip(col_names, tpl))
        json_bytes = entry_dict['json']
        stmt = Statement._from_json(json.loads(json_bytes.decode('utf-8')))
        entry_dict['mk_hash'] = stmt.get_hash()
        ret_tpl = tuple([entry_dict[col] for col in copy_col_names])
        copy_stmt_tuples.append(ret_tpl)
    return copy_stmt_tuples, copy_col_names


def _get_loaded_db(split=None, with_init_corpus=False):
    print("Creating and filling a test database:")
    db = _get_background_loaded_db()

    # Now load the statements. Much of this processing is the result of active
    # development, and once that is done, TODO: Format pickle to match
    copy_stmt_tuples, copy_col_names = _get_input_stmt_tuples()

    print("\tInserting the raw statements...")
    if split is None:
        db.copy('raw_statements', copy_stmt_tuples, copy_col_names)
        print("\tAdding agents...")
        db_util.insert_agents(db, 'raw')
        if with_init_corpus:
            print("\tAdding a preassembled corpus...")
            pa_manager = pm.PreassemblyManager()
            pa_manager.create_corpus(db)
    else:
        num_initial = int(split*len(copy_stmt_tuples))
        stmt_tuples_initial = random.sample(copy_stmt_tuples, num_initial)
        stmt_tuples_new = list(set(copy_stmt_tuples) - set(stmt_tuples_initial))
        initial_datetime = datetime.now() - timedelta(days=2)
        db.copy('raw_statements', [t + (initial_datetime,)
                                   for t in stmt_tuples_initial],
                copy_col_names + ('create_date',))
        print("\tAdding agents...")
        db_util.insert_agents(db, 'raw')
        if with_init_corpus:
            print("\tAdding a preassembled corpus from first batch of raw "
                  "stmts...")
            pa_manager = pm.PreassemblyManager(batch_size=BATCH_SIZE)
            pa_manager.create_corpus(db)
        print("\tInserting the rest of the raw statements...")
        new_datetime = datetime.now()
        db.copy('raw_statements', [t + (new_datetime,) for t in stmt_tuples_new],
                copy_col_names + ('create_date',))
        print("\tAdding agents...")
        db_util.insert_agents(db, 'raw')
    return db


def _get_stmts():
    global STMTS
    if STMTS is None:
        stmt_tuples, _ = _get_stmt_tuples()
        stmt_jsons = [json.loads(tpl[-1].decode('utf8')) for tpl in stmt_tuples]
        STMTS = stmts_from_json(stmt_jsons)
    return STMTS


def _str_large_set(s, max_num):
    if len(s) > max_num:
        values = list(s)[:max_num]
        ret_str = '{' + ', '.join([str(v) for v in values]) + ' ...}'
        ret_str += ' [length: %d]' % len(s)
    else:
        ret_str = str(s)
    return ret_str


def _do_old_fashioned_preassembly(stmts):
    grounded_stmts = ac.map_grounding(stmts)
    ms_stmts = ac.map_sequence(grounded_stmts)
    opa_stmts = ac.run_preassembly(ms_stmts, return_toplevel=False)
    return opa_stmts


def _check_against_opa_stmts(raw_stmts, pa_stmts):
    def _compare_list_elements(label, list_func, comp_func, **stmts):
        (stmt_1_name, stmt_1), (stmt_2_name, stmt_2) = list(stmts.items())
        vals_1 = [comp_func(elem) for elem in list_func(stmt_1)]
        vals_2 = []
        for element in list_func(stmt_2):
            val = comp_func(element)
            if val in vals_1:
                vals_1.remove(val)
            else:
                vals_2.append(val)
        if len(vals_1) or len(vals_2):
            print("Found mismatched %s for hash %s:\n\t%s=%s\n\t%s=%s"
                  % (label, stmt_1.get_hash(shallow=True), stmt_1_name, vals_1,
                     stmt_2_name, vals_2))
            return {'diffs': {stmt_1_name: vals_1, stmt_2_name: vals_2},
                    'stmts': {stmt_1_name: stmt_1, stmt_2_name: stmt_2}}
        return None

    opa_stmts = _do_old_fashioned_preassembly(raw_stmts)

    old_stmt_dict = {s.get_hash(shallow=True): s for s in opa_stmts}
    new_stmt_dict = {s.get_hash(shallow=True): s for s in pa_stmts}

    new_hash_set = set(new_stmt_dict.keys())
    old_hash_set = set(old_stmt_dict.keys())
    hash_diffs = {'extra_new': [new_stmt_dict[h]
                                for h in new_hash_set - old_hash_set],
                  'extra_old': [old_stmt_dict[h]
                                for h in old_hash_set - new_hash_set]}
    print(hash_diffs)
    tests = [{'funcs': {'list': lambda s: s.evidence,
                        'comp': lambda ev: ev.matches_key()},
              'label': 'evidence text',
              'results': []},
             {'funcs': {'list': lambda s: s.supports,
                        'comp': lambda s: s.get_hash(shallow=True)},
              'label': 'supports matches keys',
              'results': []},
             {'funcs': {'list': lambda s: s.supported_by,
                        'comp': lambda s: s.get_hash(shallow=True)},
              'label': 'supported-by matches keys',
              'results': []}]
    for mk_hash in (new_hash_set & old_hash_set):
        for test_dict in tests:
            res = _compare_list_elements(test_dict['label'],
                                         test_dict['funcs']['list'],
                                         test_dict['funcs']['comp'],
                                         new_stmt=new_stmt_dict[mk_hash],
                                         old_stmt=old_stmt_dict[mk_hash])
            if res is not None:
                test_dict['results'].append(res)

    # Now evaluate the results for exceptions
    assert all([len(mismatch_res) is 0
                for mismatch_res in [test_dict['results']
                                     for test_dict in tests]]),\
        ('\n'.join(['Found %d mismatches in %s.' % (len(td['results']),
                                                    td['label'])
                    for td in tests]))
    assert not any(hash_diffs.values()), "Found mismatched hashes."


@needs_py3
@attr('slow')
def test_statement_distillation():
    db = _get_loaded_db()
    assert db is not None, "Test was broken. Got None instead of db insance."
    stmts = db_util.distill_stmts(db, get_full_stmts=True)
    assert len(stmts), "Got zero statements."
    assert isinstance(list(stmts)[0], Statement), type(list(stmts)[0])
    stmt_ids = db_util.distill_stmts(db)
    assert len(stmts) == len(stmt_ids), \
        "stmts: %d, stmt_ids: %d" % (len(stmts), len(stmt_ids))
    assert isinstance(list(stmt_ids)[0], str), type(list(stmt_ids)[0])
    stmts_p = db_util.distill_stmts(db, num_procs=2)
    assert len(stmts_p) == len(stmt_ids)
    stmt_ids_p = db_util.distill_stmts(db, num_procs=2)
    assert stmt_ids_p == stmt_ids


@needs_py3
@attr('slow')
def test_preassembly_with_database():
    db = _get_loaded_db()

    # Get the set of raw statements.
    raw_stmt_list = db.select_all(db.RawStatements)
    all_raw_uuids = {raw_stmt.uuid for raw_stmt in raw_stmt_list}
    assert len(raw_stmt_list)

    # Run the preassembly initialization.
    start = datetime.now()
    pa_manager = pm.PreassemblyManager(batch_size=BATCH_SIZE)
    pa_manager.create_corpus(db)
    end = datetime.now()
    print("Duration:", end-start)
    pa_stmt_list = db.select_all(db.PAStatements)
    assert 0 < len(pa_stmt_list) < len(raw_stmt_list)
    raw_unique_link_list = db.select_all(db.RawUniqueLinks)
    assert len(raw_unique_link_list)
    all_link_uuids = {ru.raw_stmt_uuid for ru in raw_unique_link_list}
    all_link_mk_hashes = {ru.pa_stmt_mk_hash for ru in raw_unique_link_list}
    assert len(all_link_uuids - all_raw_uuids) is 0
    assert all([pa_stmt.mk_hash in all_link_mk_hashes
                for pa_stmt in pa_stmt_list])
    num_support_links = db.filter_query(db.PASupportLinks).count()
    assert num_support_links

    # Try to get all the preassembled statements from the table.
    pa_stmts = db_client.get_statements([], preassembled=True, db=db)
    assert len(pa_stmts) == len(pa_stmt_list), (len(pa_stmts),
                                                len(pa_stmt_list))

    # Now test the set of preassembled (pa) statements from the database against
    # what we get from old-fashioned preassembly (opa).
    raw_stmts = db_util.distill_stmts(db, get_full_stmts=True)
    _check_against_opa_stmts(raw_stmts, pa_stmts)


@needs_py3
@attr('slow')
def test_incremental_preassembly_with_database():
    db = _get_loaded_db(split=0.8, with_init_corpus=True)
    pa_manager = pm.PreassemblyManager(batch_size=BATCH_SIZE)
    print("Beginning supplement...")
    pa_manager.supplement_corpus(db)

    raw_stmts = db_util.distill_stmts(db, get_full_stmts=True)
    pa_stmts = db_client.get_statements([], preassembled=True, db=db)
    _check_against_opa_stmts(raw_stmts, pa_stmts)
