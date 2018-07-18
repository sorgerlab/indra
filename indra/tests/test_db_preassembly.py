from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
import random
import logging
from datetime import datetime


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
from indra.db.preassembly_manager import shash
from indra.statements import Statement
from indra.tools import assemble_corpus as ac

from nose.plugins.attrib import attr
from .util import needs_py3
from .test_db_client import _PrePaDatabaseTestSetup

from indra.statements import Phosphorylation, Agent, Evidence
from indra.db.util import NestedDict
from indra.db.util import reader_versions as rv_dict

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MAX_NUM_STMTS = 11721
BATCH_SIZE = 2017
STMTS = None

# ==============================================================================
# Support clases and functions
# ==============================================================================


def make_raw_statement_set_for_distillation():
    d = NestedDict()
    stmts = []
    target_sets = []
    bettered_sids = set()

    # Create a function which will update all possible outcome scenarios given a
    # set of some_stmts.
    def add_stmts_to_target_set(some_stmts):
        # If we don't have any target sets of statements, initialize with the
        # input statements.
        if not target_sets:
            for stmt in some_stmts:
                target_sets.append(({stmt},
                                    {stmts.index(s) for s in some_stmts
                                     if s is not stmt}))
        else:
            # Make a copy and empty the current list.
            old_target_sets = target_sets[:]
            try: # Python 3
                target_sets.clear()
            except AttributeError: # Python 2
                del target_sets[:]

            # Now for every previous scenario, pick a random possible "good"
            # statement, update the corresponding duplicate trace.
            for stmt_set, dup_set in old_target_sets:
                for stmt in some_stmts:
                    # Here we consider the possibility that each of the
                    # potential valid statements may be chosen, and record that
                    # possible alteration to the set of possible histories.
                    new_set = stmt_set.copy()
                    new_set.add(stmt)
                    new_dups = dup_set.copy()
                    new_dups |= {stmts.index(s) for s in some_stmts
                                 if s is not stmt}
                    target_sets.append((new_set, new_dups))
        return target_sets

    # Create a function to handle the creation of the metadata.
    def add_content(trid, src, tcid, reader, rv_idx, rid, a, b, ev_num, copies,
                    is_target=False):
        # Add the new statements to the over-all list.
        stmts.extend(__make_test_statements(a, b, reader, ev_num, copies))

        # If we are making multiple copies, the latest copies should have the
        # same overall hash. If it's not a copy, the hashes should be different.
        if copies > 1:
            # The above only applies if the evidence was specified to be the
            # same, otherwise it assumed the evidence, and therefore the hash,
            # is different.
            if ev_num is not None:
                assert stmts[-1].get_hash() == stmts[-2].get_hash()
            else:
                assert stmts[-1].get_hash() != stmts[-2].get_hash()

        # Populate the provenance for the dict.
        rv = rv_dict[reader][rv_idx]
        r_dict = d[trid][src][tcid][reader][rv][rid]

        # If the evidence variation was specified, the evidence in any copies is
        # identical, and they will all have the same hash. Else, the hash is
        # different and the statements need to be iterated over.
        if ev_num is not None:
            s_hash = stmts[-1].get_hash()

            # Check to see if we have a matching statment yet.
            if r_dict.get(s_hash) is None:
                r_dict[s_hash] = set()

            # Set the value
            d[trid][src][tcid][reader][rv][rid][stmts[-1].get_hash()] |= \
                {(stmts.index(s), s) for s in stmts[-copies:]}
        else:
            for s in stmts[-copies:]:
                s_hash = s.get_hash()
                if r_dict.get(s_hash) is None:
                    r_dict[s_hash] = set()
                d[trid][src][tcid][reader][rv][rid][s_hash].add(
                    (stmts.index(s), s)
                    )

        # If this/these statement/s is intended to be picked up, add it/them to
        # the target sets.
        if is_target:
            global target_sets
            target_sets = add_stmts_to_target_set(stmts[-copies:])
        return

    # We produced statements a coupld of times with and old reader version
    #           trid         tcid        reader vrsn idx   distinct evidence id
    #           |  source    |  reader   |  reading id     |  number of copies
    #           |  |         |  |        |  |  Agents      |  |  Is it a target?
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A1', 'B1', 1, 2, False)
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A1', 'B1', 2, 1)
    add_content(1, 'pubmed', 1, 'reach', 0, 1, 'A2', 'B2', 1, 1)

    # Do it again for a new reader version.
    add_content(1, 'pubmed', 1, 'reach', 1, 2, 'A1', 'B1', 1, 2)
    add_content(1, 'pubmed', 1, 'reach', 1, 2, 'A1', 'B1', 2, 1)

    # Add some for sparser.
    add_content(1, 'pubmed', 1, 'sparser', 1, 3, 'A1', 'B1', 1, 2)
    add_content(1, 'pubmed', 1, 'sparser', 1, 3, 'A2', 'B2', 1, 1)

    # Now add statements from another source.
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A1', 'B1', 1, 2)
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A1', 'B1', 2, 1)
    add_content(1, 'pmc_oa', 2, 'reach', 0, 4, 'A2', 'B2', 1, 1)

    # All the statements up until now will be skipped, if all goes well.
    bettered_sids |= set(range(len(stmts)))

    # ...and again for a new reader version.
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A1', 'B1', 1, 2, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A1', 'B1', 2, 1, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A2', 'B2', 1, 1, True)
    add_content(1, 'pmc_oa', 2, 'reach', 1, 4, 'A3', 'B3', 1, 1, True)

    # Add some results from sparser
    add_content(1, 'pmc_oa', 2, 'sparser', 1, 5, 'A1', 'B1', 1, 2, True)
    add_content(1, 'pmc_oa', 2, 'sparser', 1, 5, 'A2', 'B2', 1, 1, True)

    # Add some content for another text ref.
    add_content(2, 'pmc_oa', 3, 'sparser', 1, 6, 'A3', 'B3', 1, 1, True)
    add_content(2, 'manuscripts', 4, 'sparser', 1, 7, 'A3', 'B3', 1, 1)

    # This last statement should also be skipped, if all goes well.
    bettered_sids.add(len(stmts) - 1)

    return d, stmts, target_sets, bettered_sids


def __make_test_statements(a, b, source_api, ev_num=None, copies=1):
    stmts = []
    A = Agent(a)
    B = Agent(b)
    for i in range(copies):
        if ev_num is None:
            ev_num = i
        ev_text = "Evidence %d for %s phosphorylates %s." % (ev_num, a, b)
        ev_list = [Evidence(text=ev_text, source_api=source_api)]
        stmts.append(Phosphorylation(Agent(A), Agent(B), evidence=ev_list))
    return stmts


class _DatabaseTestSetup(_PrePaDatabaseTestSetup):
    """This object is used to setup the test database into various configs."""
    def add_statements(self, fraction=1, pam=None):
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
            num_stmts = int(fraction*len(available_tuples))
            input_tuples = random.sample(available_tuples, num_stmts)
        else:
            input_tuples = available_tuples

        self.insert_the_statements(input_tuples)

        if pam:
            print("Preassembling new statements...")
            if self.used_stmt_tuples:
                pam.supplement_corpus(self.test_db)
            else:
                pam.create_corpus(self.test_db)

        self.used_stmt_tuples |= set(input_tuples)
        return


def _get_loaded_db(num_stmts, split=None, pam=None):
    print("Creating and filling a test database:")
    dts = _DatabaseTestSetup(num_stmts)
    dts.load_background()

    if split is None:
        dts.add_statements(pam=pam)
    else:
        dts.add_statements(split, pam=pam)
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


def _get_opa_input_stmts(db):
    stmt_nd = db_util._get_reading_statement_dict(db, get_full_stmts=True)
    reading_stmts, _, _ =\
        db_util._get_filtered_rdg_statements(stmt_nd, get_full_stmts=True,
                                             ignore_duplicates=True)
    db_stmts = db_client.get_statements([db.RawStatements.reading_id == None],
                                        preassembled=False, db=db)
    stmts = reading_stmts | set(db_stmts)
    print("Got %d statements for opa." % len(stmts))
    return stmts


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
                  % (label, shash(stmt_1), stmt_1_name, vals_1, stmt_2_name,
                     vals_2))
            return {'diffs': {stmt_1_name: vals_1, stmt_2_name: vals_2},
                    'stmts': {stmt_1_name: stmt_1, stmt_2_name: stmt_2}}
        return None

    opa_stmts = _do_old_fashioned_preassembly(raw_stmts)

    old_stmt_dict = {shash(s): s for s in opa_stmts}
    new_stmt_dict = {shash(s): s for s in pa_stmts}

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
    tests = [{'funcs': {'list': lambda s: s.evidence[:],
                        'comp': lambda ev: ev.matches_key()},
              'label': 'evidence text',
              'results': []},
             {'funcs': {'list': lambda s: s.supports[:],
                        'comp': lambda s: shash(s)},
              'label': 'supports matches keys',
              'results': []},
             {'funcs': {'list': lambda s: s.supported_by[:],
                        'comp': lambda s: shash(s)},
              'label': 'supported-by matches keys',
              'results': []}]
    comp_hashes = new_hash_set & old_hash_set
    for mk_hash in comp_hashes:
        for test_dict in tests:
            res = _compare_list_elements(test_dict['label'],
                                         test_dict['funcs']['list'],
                                         test_dict['funcs']['comp'],
                                         new_stmt=new_stmt_dict[mk_hash],
                                         old_stmt=old_stmt_dict[mk_hash])
            if res is not None:
                test_dict['results'].append(res)

    def all_tests_passed():
        test_results = [not any(hash_diffs.values())]
        for td in tests:
            test_results.append(len(td['results']) == 0)
        print("%d/%d tests passed." % (sum(test_results), len(test_results)))
        return all(test_results)

    def write_report(num_comps):
        ret_str = "Some tests failed:\n"
        ret_str += ('Found %d/%d extra old stmts and %d/%d extra new stmts.\n'
                    % (len(hash_diffs['extra_old']), len(old_hash_set),
                       len(hash_diffs['extra_new']), len(new_hash_set)))
        for td in tests:
            ret_str += ('Found %d/%d mismatches in %s.\n'
                        % (len(td['results']), num_comps, td['label']))
        return ret_str

    # Now evaluate the results for exceptions
    assert all_tests_passed(), write_report(len(comp_hashes))


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
                 % (str(s), o.id, o.uuid[:8] + '...', o.type,
                    o.indra_version[:14] + '...', o.mk_hash))
        if other_stmt_keys and shash(s) in other_stmt_keys:
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
              % (s.uuid, s.get_hash(), shash(s)))
        print('-'*100)
        db_pas = db.select_one(db.PAStatements,
                               db.PAStatements.mk_hash == shash(s))
        print('\tPA statement:', db_pas.__dict__ if db_pas else '~')
        print('-'*100)
        db_s = db.select_one(db.RawStatements, db.RawStatements.uuid == s.uuid)
        print('\tRaw statement:', str_imp(db_s, uuid, other_stmt_keys))
        if db_s is None:
            continue
        print('-'*100)
        if db_s.reading_id is None:
            print("Statement was from a database: %s" % db_s.db_info_id)
            continue
        db_r = db.select_one(db.Reading, db.Reading.id == db_s.reading_id)
        print('\tReading:', str_imp(db_r))
        tc = db.select_one(db.TextContent,
                           db.TextContent.id == db_r.text_content_id)
        print('\tText Content:', str_imp(tc))
        tr = db.select_one(db.TextRef, db.TextRef.id == tc.text_ref_id)
        print('\tText ref:', str_imp(tr))
        print('-'*100)
        for tc in db.select_all(db.TextContent,
                                db.TextContent.text_ref_id == tr.id):
            print('\t', str_imp(tc))
            for r in db.select_all(db.Reading,
                                   db.Reading.text_content_id == tc.id):
                print('\t\t', str_imp(r))
                for s in db.select_all(db.RawStatements,
                                       db.RawStatements.reading_id == r.id):
                    print('\t\t\t', str_imp(s, uuid, other_stmt_keys))
        print('='*100)


# ==============================================================================
# Generic test definitions
# ==============================================================================


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


@needs_py3
def _check_preassembly_with_database(num_stmts, batch_size, n_proc=1):
    db = _get_loaded_db(num_stmts)

    # Now test the set of preassembled (pa) statements from the database against
    # what we get from old-fashioned preassembly (opa).
    opa_inp_stmts = _get_opa_input_stmts(db)

    # Get the set of raw statements.
    raw_stmt_list = db.select_all(db.RawStatements)
    all_raw_ids = {raw_stmt.id for raw_stmt in raw_stmt_list}
    assert len(raw_stmt_list)

    # Run the preassembly initialization.
    start = datetime.now()
    pa_manager = pm.PreassemblyManager(batch_size=batch_size, n_proc=n_proc,
                                       print_logs=True)
    pa_manager.create_corpus(db)
    end = datetime.now()
    print("Duration:", end-start)

    # Make sure the number of pa statements is within reasonable bounds.
    pa_stmt_list = db.select_all(db.PAStatements)
    assert 0 < len(pa_stmt_list) < len(raw_stmt_list)

    # Check the evidence links.
    raw_unique_link_list = db.select_all(db.RawUniqueLinks)
    assert len(raw_unique_link_list)
    all_link_ids = {ru.raw_stmt_id for ru in raw_unique_link_list}
    all_link_mk_hashes = {ru.pa_stmt_mk_hash for ru in raw_unique_link_list}
    assert len(all_link_ids - all_raw_ids) is 0
    assert all([pa_stmt.mk_hash in all_link_mk_hashes
                for pa_stmt in pa_stmt_list])

    # Check the support links.
    sup_links = db.select_all([db.PASupportLinks.supporting_mk_hash,
                               db.PASupportLinks.supported_mk_hash])
    assert sup_links
    assert not any([l[0] == l[1] for l in sup_links]),\
        "Found self-support in the database."

    # Try to get all the preassembled statements from the table.
    pa_stmts = db_client.get_statements([], preassembled=True, db=db,
                                        with_support=True)
    assert len(pa_stmts) == len(pa_stmt_list), (len(pa_stmts),
                                                len(pa_stmt_list))

    self_supports = {
        shash(s): shash(s) in {shash(s_) for s_ in s.supported_by + s.supports}
        for s in pa_stmts
        }
    if any(self_supports.values()):
        assert False, "Found self-support in constructed pa statement objects."

    _check_against_opa_stmts(db, opa_inp_stmts, pa_stmts)
    return


@needs_py3
def _check_db_pa_supplement(num_stmts, batch_size, split=0.8, n_proc=1):
    pa_manager = pm.PreassemblyManager(batch_size=batch_size, n_proc=n_proc,
                                       print_logs=True)
    db = _get_loaded_db(num_stmts, split=split, pam=pa_manager)
    opa_inp_stmts = _get_opa_input_stmts(db)
    start = datetime.now()
    print("Beginning supplement...")
    pa_manager.supplement_corpus(db)
    end = datetime.now()
    print("Duration of incremental update:", end-start)

    pa_stmts = db_client.get_statements([], preassembled=True, db=db,
                                        with_support=True)
    _check_against_opa_stmts(db, opa_inp_stmts, pa_stmts)
    return


# ==============================================================================
# Specific Tests
# ==============================================================================


def test_distillation_on_curated_set():
    stmt_dict, stmt_list, target_sets, target_bettered_ids = \
        make_raw_statement_set_for_distillation()
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
        make_raw_statement_set_for_distillation()
    filtered_id_set, duplicate_ids, bettered_ids = \
        db_util._get_filtered_rdg_statements(stmt_dict, get_full_stmts=False)
    assert len(filtered_id_set) == len(filtered_set), \
        (len(filtered_set), len(filtered_id_set))


@attr('nonpublic')
def test_statement_distillation_small():
    _check_statement_distillation(1000)


@attr('nonpublic', 'slow')
def test_statement_distillation_large():
    _check_statement_distillation(11721)


@attr('nonpublic', 'slow')
def test_statement_distillation_extra_large():
    _check_statement_distillation(1001721)


@attr('nonpublic')
def test_db_preassembly_small():
    _check_preassembly_with_database(200, 40)


@attr('nonpublic', 'slow')
def test_db_preassembly_large():
    _check_preassembly_with_database(11721, 2017)


@attr('nonpublic', 'slow')
def test_db_preassembly_extra_large():
    _check_preassembly_with_database(101721, 20017)


@attr('nonpublic', 'slow')
def test_db_preassembly_supremely_large():
    _check_preassembly_with_database(1001721, 200017)


@attr('nonpublic')
def test_db_incremental_preassembly_small():
    _check_db_pa_supplement(200, 40)


@attr('nonpublic', 'slow')
def test_db_incremental_preassembly_large():
    _check_db_pa_supplement(11721, 2017)


@attr('nonpublic', 'slow')
def test_db_incremental_preassembly_very_large():
    _check_db_pa_supplement(100000, 20000, n_proc=2)


@attr('nonpublic', 'slow')
def test_db_incremental_preassembly_1M():
    _check_db_pa_supplement(1000000, 200000, n_proc=6)
