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
from indra.db import preassembly_script as pas
from indra.statements import stmts_from_json, Statement
from indra.tools import assemble_corpus as ac

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
STMT_PICKLE_FILE = os.path.join(THIS_DIR, 'test_stmts_tuples.pkl')
MAX_NUM_STMTS = 10000
STMTS = None


def _load_tuples(fname):
    with open(os.path.join(THIS_DIR, fname), 'rb') as f:
        ret_tuples = pickle.load(f)
    return ret_tuples


def _get_stmt_tuples():
    with open(STMT_PICKLE_FILE, 'rb') as f:
        stmt_tuples = pickle.load(f)
    col_names = stmt_tuples.pop(0)
    if len(stmt_tuples) > MAX_NUM_STMTS:
        stmt_tuples = random.sample(stmt_tuples, MAX_NUM_STMTS)
    return stmt_tuples, col_names


def _get_loaded_db(split=None, with_init_corpus=False):
    print("Creating and filling a test database:")
    db = db_util.get_test_db()
    db._clear(force=True)
    stmt_tuples, col_names = _get_stmt_tuples()

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

    # Now load the statements. Much of this processing is the result of active
    # development, and once that is done, TODO: Format pickle to match
    print("\tPrepping the raw statements...")
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

    print("\tInserting the raw statements...")
    if split is None:
        db.copy('raw_statements', copy_stmt_tuples, copy_col_names)
        if with_init_corpus:
            print("\tAdding a preassembled corpus...")
            pa_manager = pas.PreassemblyManager()
            pa_manager.create_corpus(db)
    else:
        num_initial = int(split*len(copy_stmt_tuples))
        stmt_tuples_initial = random.sample(copy_stmt_tuples, num_initial)
        stmt_tuples_new = list(set(copy_stmt_tuples) - set(stmt_tuples_initial))
        initial_datetime = datetime.now() - timedelta(days=2)
        db.copy('raw_statements', [t + (initial_datetime,)
                                   for t in stmt_tuples_initial],
                copy_col_names + ('create_date',))
        if with_init_corpus:
            print("\tAdding a preassembled corpus from first batch of raw "
                  "stmts...")
            pa_manager = pas.PreassemblyManager()
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
    opa_stmts = _do_old_fashioned_preassembly(raw_stmts)

    old_stmt_dict = {s.get_shallow_hash(): s for s in opa_stmts}
    new_stmt_dict = {s.get_shallow_hash(): s for s in pa_stmts}

    new_hash_set = set(new_stmt_dict.keys())
    old_hash_set = set(old_stmt_dict.keys())
    assert new_hash_set == old_hash_set, \
        (new_hash_set - old_hash_set, old_hash_set - new_hash_set)
    mismatched_evidence_texts = []
    for mk_hash in new_hash_set:
        new_ev_texts = [ev.text for ev in new_stmt_dict[mk_hash].evidence]
        old_ev_texts = []
        for ev in old_stmt_dict[mk_hash].evidence:
            if ev.text in new_ev_texts:
                new_ev_texts.remove(ev.text)
            else:
                old_ev_texts.append(ev.text)
        if len(old_ev_texts) or len(new_ev_texts):
            print("Found mismatched evidence: old: %s, new: %s."
                  % (old_ev_texts, new_ev_texts))
            mismatched_evidence_texts.append({'old': old_ev_texts,
                                              'new': new_ev_texts})
    assert not len(mismatched_evidence_texts), \
        "Found %d mismatched evidence." % len(mismatched_evidence_texts)


def test_preassembly_without_database():
    stmts = _get_stmts()
    unique_stmt_dict, evidence_links, support_links = \
        pas.process_statements(stmts)
    assert len(unique_stmt_dict)
    total_evidence = len(evidence_links)
    assert len(unique_stmt_dict) <= total_evidence <= len(stmts), \
        ("Got %d ev links for %d stmts and %d unique statements (should be "
         "between)." % (total_evidence, len(stmts), len(unique_stmt_dict)))
    num_unique_stmt_keys = len({t[0] for t in evidence_links})
    assert num_unique_stmt_keys == len(unique_stmt_dict), \
        ("Got %d ev sets for %d unique stmts."
         % (num_unique_stmt_keys, len(unique_stmt_dict)))
    assert len(support_links)

    opa_stmts = _do_old_fashioned_preassembly(stmts)
    new_hash_set = set(unique_stmt_dict.keys())
    old_hash_set = {s.get_shallow_hash() for s in opa_stmts}
    assert new_hash_set == old_hash_set, \
        (new_hash_set - old_hash_set, old_hash_set - new_hash_set)
    return


def test_incremental_preassmbly_without_database():
    stmts = _get_stmts()

    # For comparison, preassemble the entire corpus.
    full_unique_stmts, full_evidence_links, full_support_links = \
        pas.process_statements(stmts)

    # Randomly split the sample 80-20
    init_stmts = random.sample(stmts, int(0.8*len(stmts)))
    new_stmts = list(set(stmts) - set(init_stmts))

    # Run preassmbly on the "init" corpus (the 80)
    init_unique_stmts, init_evidence_links, init_support_links = \
        pas.process_statements(init_stmts)
    assert len(init_support_links)

    # Make sure the "new" statements actually have at least some links to add
    _, new_only_ev_links, new_only_mk_links = pas.process_statements(new_stmts)
    assert len(new_only_mk_links), "Test not useful without new matches."
    print("Evidence links from new stmts:", len(new_only_ev_links))
    print("Support links from new stmts:", len(new_only_mk_links))

    # Now merge in the "new" statements (the 20)
    new_unique_stmt_dict, new_evidence_links, new_support_links = \
        pas.get_increment_links(init_unique_stmts, new_stmts)
    updated_unique_stmts, updated_evidence_links, updated_support_links = \
        pas.merge_statements(init_unique_stmts, init_evidence_links,
                             init_support_links, new_unique_stmt_dict,
                             new_evidence_links, new_support_links)

    # Check that we got all the same statements (trivial)
    assert len(updated_unique_stmts) == len(full_unique_stmts), \
        ("Got %d unique stmts from update, but %d from pre-assembly of all "
         "stmts." % (len(updated_unique_stmts), len(full_unique_stmts)))

    # Check that the evidence matches up (easy)
    # fevl_set = pas.flatten_evidence_dict(full_evidence_links)
    # uevl_set = pas.flatten_evidence_dict(updated_evidence_links)
    missed_evidence_links = full_evidence_links - updated_evidence_links
    extra_evidence_links = updated_evidence_links - full_evidence_links
    assert not len(missed_evidence_links) and not len(extra_evidence_links), \
        ("Some evidence links missed: %s, and/or some evidence links added: %s"
         % (_str_large_set(missed_evidence_links, 5),
            _str_large_set(extra_evidence_links, 5)))

    # Check that we have the same match key links (less easy)
    missed_supports = full_support_links - updated_support_links
    extra_supports = updated_support_links - full_support_links
    assert not len(extra_supports) and not len(missed_supports), \
        ("Some match key links missed: %s and/or some match key links added: %s"
         % (missed_supports, extra_supports))
    return


def test_statement_distillation():
    db = _get_loaded_db()
    assert db is not None, "Test was broken. Got None instead of db insance."
    stmt_nd, stmts = db_util.distill_stmts_from_reading(db, get_full_stmts=True)
    assert len(stmt_nd.keys()), "Got no nested dict."
    assert len(stmts), "Got zero statements."
    assert isinstance(list(stmts)[0], Statement), type(list(stmts)[0])
    stmt_id_nd, stmt_ids = db_util.distill_stmts_from_reading(db)
    assert len(stmt_id_nd) == len(stmt_nd), \
        "stmt_id_nd: %d, stmt_nd: %d" % (len(stmt_id_nd), len(stmt_nd))
    assert len(stmts) == len(stmt_ids), \
        "stmts: %d, stmt_ids: %d" % (len(stmts), len(stmt_ids))
    assert isinstance(list(stmt_ids)[0], str), type(list(stmt_ids)[0])


def test_preassembly_with_database():
    db = _get_loaded_db()

    # Get the set of raw statements.
    raw_stmt_list = db.select_all(db.RawStatements)
    all_raw_uuids = {raw_stmt.uuid for raw_stmt in raw_stmt_list}
    assert len(raw_stmt_list)

    # Run the preassembly initialization.
    pa_manager = pas.PreassemblyManager()
    pa_manager.create_corpus(db)
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
    pa_stmts = db_util.get_statements([], preassembled=True, db=db)
    assert len(pa_stmts) == len(pa_stmt_list), (len(pa_stmts),
                                                len(pa_stmt_list))

    # Now test the set of preassembled (pa) statements from the database against
    # what we get from old-fashioned preassembly (opa).
    raw_stmts = db_util.make_stmts_from_db_list(raw_stmt_list)
    _check_against_opa_stmts(raw_stmts, pa_stmts)


def test_incremental_preassembly_with_database():
    db = _get_loaded_db(split=0.8, with_init_corpus=True)
    pa_manager = pas.PreassemblyManager()
    print("Beginning supplement...")
    pa_manager.supplement_corpus(db)

    raw_stmts = db_util.get_statements([], preassembled=False, db=db)
    pa_stmts = db_util.get_statements([], preassembled=True, db=db)
    _check_against_opa_stmts(raw_stmts, pa_stmts)
