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
from .make_raw_statement_test_set import make_raw_statement_test_set

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MAX_NUM_STMTS = 11721
BATCH_SIZE = 2017
STMTS = None


class _DatabaseTestSetup(object):
    """This object is used to setup the test database into various configs."""
    def __init__(self, max_total_stmts):
        self.test_db = db_util.get_test_db()
        self.test_db._clear(force=True)
        with open(os.path.join(THIS_DIR, 'db_pa_test_input_1M.pkl'), 'rb') as f:
            self.test_data = pickle.load(f)

        if max_total_stmts < len(self.test_data['raw_statements']['tuples']):
            self.stmt_tuples = random.sample(
                self.test_data['raw_statements']['tuples'],
                max_total_stmts
                )
        else:
            self.stmt_tuples = self.test_data['raw_statements']['tuples']

        self.used_stmt_tuples = set()
        return

    def get_available_stmt_tuples(self):
        return list(set(self.stmt_tuples) - self.used_stmt_tuples)

    def load_background(self):
        """Load in all the background provenance metadata (e.g. text_ref).

        Note: This must be done before you try to load any statements.
        """
        for tbl in ['text_ref', 'text_content', 'reading', 'db_info']:
            print("Loading %s..." % tbl)
            self.test_db.copy(tbl, self.test_data[tbl]['tuples'],
                              self.test_data[tbl]['cols'])
        return

    def add_statements(self, fraction=1, with_pa=False):
        """Add statements and agents to the database.

        Parameters
        ----------
        fraction : float between 0 and 1
            Default is 1. The fraction of remaining statements to be added.
        with_pa : bool
            Default False. Choose to run pre-assembly/incremental-preassembly
            on the added statements.
        """
        available_tuples = self.get_available_stmt_tuples()
        if fraction is not 1:
            num_stmts = fraction*len(available_tuples)
            input_tuples = random.sample(available_tuples, num_stmts)
        else:
            input_tuples = available_tuples

        print("Loading %d statements..." % len(input_tuples))
        if hasattr(self.test_db.RawStatements, 'id'):
            self.test_db.copy('raw_statements', input_tuples,
                               self.test_data['raw_statements']['cols'])
        else:
            self.test_db.copy('raw_statements', [t[1:] for t in input_tuples],
                              self.test_data['raw_statements']['cols'][1:])

        print("Inserting agents...")
        db_util.insert_agents(self.test_db, 'raw')

        if with_pa:
            print("Preassembling new statements...")
            if len(input_tuples) > 100:
                batch_size = len(input_tuples)//10
                pam = pm.PreassemblyManager(1, batch_size)
            else:
                pam = pm.PreassemblyManager()

            if self.used_stmt_tuples:
                pam.supplement_corpus(self.test_db)
            else:
                pam.create_corpus(self.test_db)

        return


def _get_loaded_db(num_stmts, split=None, with_init_corpus=False):
    print("Creating and filling a test database:")
    dts = _DatabaseTestSetup(num_stmts)
    dts.load_background()

    if split is None:
        dts.add_statements(with_pa=with_init_corpus)
    else:
        dts.add_statements(split, with_pa=with_init_corpus)
        dts.add_statements()
    return dts.test_db


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


def _check_against_opa_stmts(db, raw_stmts, pa_stmts):
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
    if hash_diffs['extra_new']:
        elaborate_on_hash_diffs(db, 'new', hash_diffs['extra_new'],
                                old_stmt_dict.keys())
    if hash_diffs['extra_old']:
        elaborate_on_hash_diffs(db, 'old', hash_diffs['extra_old'],
                                new_stmt_dict.keys())
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


def str_imp(o, uuid=None, other_stmt_keys=None):
    if o is None:
        return '~'
    cname = o.__class__.__name__
    if cname == 'TextRef':
        return ('<TextRef: trid: %s, pmid: %s, pmcid: %s>'
                % (o.id, o.pmid, o.pmcid))
    if cname == 'TextContent':
        return ('<TextContent: tcid: %s, trid: %s, src: %s>'
                % (o.id, o.text_ref_id, o.source))
    if cname == 'Reading':
        return ('<Reading: rid: %s, tcid: %s, reader: %s, rv: %s>'
                % (o.id, o.text_content_id, o.reader, o.reader_version))
    if cname == 'RawStatements':
        s = Statement._from_json(json.loads(o.json.decode()))
        s_str = ('<RawStmt: %s sid: %s, uuid: %s, type: %s, iv: %s, hash: %s>'
                 % (str(s), o.id, o.uuid, o.type, o.indra_version, o.mk_hash))
        if other_stmt_keys and s.get_hash(shallow=True) in other_stmt_keys:
            s_str = '+' + s_str
        if s.uuid == uuid:
            s_str = '*' + s_str
        return s_str


def elaborate_on_hash_diffs(db, lbl, stmt_list, other_stmt_keys):
    print("#"*100)
    print("Elaboration on extra %s statements:" % lbl)
    print("#"*100)
    for s in stmt_list:
        print(s)
        uuid = s.uuid
        print('-'*100)
        print('uuid: %s\nhash: %s\nshallow hash: %s'
              % (s.uuid, s.get_hash(), s.get_hash(shallow=True)))
        print('-'*100)
        db_pas = db.select_one(db.PAStatements,
                               db.PAStatements.mk_hash == s.get_hash(shallow=True))
        print('\tPA statement:', db_pas.__dict__ if db_pas else '~')
        print('-'*100)
        db_s = db.select_one(db.RawStatements, db.RawStatements.uuid == s.uuid)
        print('\tRaw statement:', str_imp(db_s, uuid, other_stmt_keys))
        if db_s is None:
            continue
        print('-'*100)
        db_r = db.select_one(db.Reading, db.Reading.id == db_s.reading_id)
        print('\tReading:', str_imp(db_r))
        print('\tOther Raw Statements:')
        for s in db.select_all(db.RawStatements,
                               db.RawStatements.reading_id == db_r.id):
            print('\t\t', str_imp(s, uuid, other_stmt_keys))
        print('-'*100)
        tc = db.select_one(db.TextContent,
                           db.TextContent.id == db_r.text_content_id)
        print('\tText Content:', str_imp(tc))
        print('\tOther Readings:')
        for r in db.select_all(db.Reading, db.Reading.text_content_id == tc.id):
            print('\t\t', str_imp(r))
            for s in db.select_all(db.RawStatements,
                                   db.RawStatements.reading_id == r.id):
                print('\t\t\t', str_imp(s, uuid, other_stmt_keys))
        print('-'*100)
        tr = db.select_one(db.TextRef, db.TextRef.id == tc.text_ref_id)
        print('\tText ref:', str_imp(tr))
        print('\tOther Content:')
        for tc in db.select_all(db.TextContent,
                                db.TextContent.text_ref_id == tr.id):
            print('\t\t ', str_imp(tc))
            for r in db.select_all(db.Reading,
                                   db.Reading.text_content_id == tc.id):
                print('\t\t\t', str_imp(r))
                for s in db.select_all(db.RawStatements,
                                       db.RawStatements.reading_id == r.id):
                    print('\t\t\t\t', str_imp(s, uuid, other_stmt_keys))
        print('='*100)


def test_distillation_on_curated_set():
    stmt_dict, stmt_list, target_sets, target_bettered_ids = \
        make_raw_statement_test_set()
    filtered_set, duplicate_ids, bettered_ids = \
        db_util._get_filtered_rdg_statements(stmt_dict, get_full_stmts=True)
    for stmt_set, dup_set in target_sets:
        if stmt_set == filtered_set:
            break
    else:
        assert False, "Filtered set does not match any valid possibilities."
    assert bettered_ids == target_bettered_ids
    assert dup_set == duplicate_ids, (dup_set - duplicate_ids,
                                      duplicate_ids - dup_set)
    stmt_dict, stmt_list, target_sets, target_bettered_ids = \
        make_raw_statement_test_set()
    filtered_id_set, duplicate_ids, bettered_ids = \
        db_util._get_filtered_rdg_statements(stmt_dict, get_full_stmts=False)
    assert len(filtered_id_set) == len(filtered_set), \
        (len(filtered_set), len(filtered_id_set))


@needs_py3
def _check_statement_distillation(num_stmts):
    db = _get_loaded_db(num_stmts)
    assert db is not None, "Test was broken. Got None instead of db insance."
    stmts = db_util.distill_stmts(db, get_full_stmts=True)
    assert len(stmts), "Got zero statements."
    assert isinstance(list(stmts)[0], Statement), type(list(stmts)[0])
    stmt_ids = db_util.distill_stmts(db)
    assert len(stmts) == len(stmt_ids), \
        "stmts: %d, stmt_ids: %d" % (len(stmts), len(stmt_ids))
    assert isinstance(list(stmt_ids)[0], int), type(list(stmt_ids)[0])
    stmts_p = db_util.distill_stmts(db, num_procs=2)
    assert len(stmts_p) == len(stmt_ids)
    stmt_ids_p = db_util.distill_stmts(db, num_procs=2)
    assert stmt_ids_p == stmt_ids


@attr('nonpublic')
def test_statement_distillation_small():
    _check_statement_distillation(1000)


@attr('nonpublic', 'slow')
def test_statement_distillation_large():
    _check_statement_distillation(11721)


@attr('nonpublic', 'slow')
def test_statement_distillation_extra_large():
    _check_statement_distillation(1001721)


@needs_py3
def _check_preassembly_with_database(num_stmts, batch_size):
    db = _get_loaded_db(num_stmts)

    # Get the set of raw statements.
    raw_stmt_list = db.select_all(db.RawStatements)
    all_raw_ids = {raw_stmt.id for raw_stmt in raw_stmt_list}
    assert len(raw_stmt_list)

    # Run the preassembly initialization.
    start = datetime.now()
    pa_manager = pm.PreassemblyManager(batch_size=batch_size)
    pa_manager.create_corpus(db)
    end = datetime.now()
    print("Duration:", end-start)
    pa_stmt_list = db.select_all(db.PAStatements)
    assert 0 < len(pa_stmt_list) < len(raw_stmt_list)
    raw_unique_link_list = db.select_all(db.RawUniqueLinks)
    assert len(raw_unique_link_list)
    all_link_ids = {ru.raw_stmt_id for ru in raw_unique_link_list}
    all_link_mk_hashes = {ru.pa_stmt_mk_hash for ru in raw_unique_link_list}
    assert len(all_link_ids - all_raw_ids) is 0
    assert all([pa_stmt.mk_hash in all_link_mk_hashes
                for pa_stmt in pa_stmt_list])
    num_support_links = db.filter_query(db.PASupportLinks).count()
    assert num_support_links

    # Try to get all the preassembled statements from the table.
    pa_stmts = db_client.get_statements([], preassembled=True, db=db,
                                        with_support=True)
    assert len(pa_stmts) == len(pa_stmt_list), (len(pa_stmts),
                                                len(pa_stmt_list))

    # Now test the set of preassembled (pa) statements from the database against
    # what we get from old-fashioned preassembly (opa).
    raw_stmts = db_client.get_statements([], preassembled=False, db=db)
    _check_against_opa_stmts(db, raw_stmts, pa_stmts)


@attr('nonpublic')
def test_db_preassembly_small():
    _check_preassembly_with_database(200, 40)


@attr('nonpublic', 'slow')
def test_db_preassembly_large():
    _check_preassembly_with_database(11721, 2017)


@attr('nonpublic', 'slow')
def test_db_preassembly_extra_large():
    _check_preassembly_with_database(101721, 2017)


@attr('nonpublic', 'slow')
def test_db_preassembly_supremely_large():
    _check_preassembly_with_database(1001721, 200017)

@needs_py3
def _check_db_pa_supplement(num_stmts, batch_size, split=0.8):
    db = _get_loaded_db(num_stmts, split=split, with_init_corpus=True)
    start = datetime.now()
    pa_manager = pm.PreassemblyManager(batch_size=batch_size)
    print("Beginning supplement...")
    pa_manager.supplement_corpus(db)
    end = datetime.now()
    print("Duration of incremental update:", end-start)

    raw_stmts = db_client.get_statements([], preassembled=False, db=db)
    pa_stmts = db_client.get_statements([], preassembled=True, db=db,
                                        with_support=True)
    _check_against_opa_stmts(db, raw_stmts, pa_stmts)


@attr('nonpublic')
def test_db_incremental_preassembly_small():
    _check_db_pa_supplement(200, 40)


@attr('nonpublic', 'slow')
def test_db_incremental_preassembly_large():
    _check_db_pa_supplement(11721, 2017)


@attr('nonpublic', 'slow')
def test_db_incremental_preassembly_very_large():
    _check_db_pa_supplement(100000, 20000)
