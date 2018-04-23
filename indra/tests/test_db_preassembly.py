from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
import pickle
import random
import logging

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
from indra.statements import stmts_from_json


STMT_PICKLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'test_stmts_tuples.pkl')
MAX_NUM_STMTS = 10000
STMTS = None


def _get_stmts():
    global STMTS
    if STMTS is None:
        with open(STMT_PICKLE_FILE, 'rb') as f:
            stmt_tuples = pickle.load(f)
        col_names = stmt_tuples.pop(0)
        stmt_jsons = [json.loads(tpl[-1].decode('utf8')) for tpl in stmt_tuples]
        if len(stmt_jsons) <= MAX_NUM_STMTS:
            STMTS = stmts_from_json(stmt_jsons)
        else:
            STMTS = stmts_from_json(random.sample(stmt_jsons, MAX_NUM_STMTS))
    return STMTS


def _str_large_set(s, max_num):
    if len(s) > max_num:
        values = list(s)[:max_num]
        ret_str = '{' + ', '.join([str(v) for v in values]) + ' ...}'
        ret_str += ' [length: %d]' % len(s)
    else:
        ret_str = str(s)
    return ret_str


def test_preassembly_without_database():
    stmts = _get_stmts()
    unique_stmt_dict, evidence_links, match_key_links = \
        pas.process_statements(stmts)
    assert len(unique_stmt_dict)
    total_evidence = len(pas.flatten_evidence_dict(evidence_links))
    assert total_evidence == len(stmts), \
        "Got %d ev links for %d stmts." % (total_evidence, len(stmts))
    assert len(evidence_links) == len(unique_stmt_dict), \
        ("Got %d ev sets for %d unique stmts."
         % (len(evidence_links), len(unique_stmt_dict)))
    assert len(match_key_links)


def test_incremental_preassmbly_without_database():
    stmts = _get_stmts()

    # For comparison, preassemble the entire corpus.
    full_unique_stmts, full_evidence_links, full_match_key_links = \
        pas.process_statements(stmts)

    # Randomly split the sample 80-20
    original_stmts = random.sample(stmts, int(0.8*len(stmts)))
    new_stmts = list(set(stmts) - set(original_stmts))

    # Run preassmbly on the "original" corpus (the 80)
    original_unique_stmts, orignal_evidence_links, original_match_key_links = \
        pas.process_statements(original_stmts)
    assert len(original_match_key_links)

    # Make sure the "new" statements actually have at least some links to add
    _, new_only_ev_links, new_only_mk_links = pas.process_statements(new_stmts)
    assert len(new_only_mk_links), "Test not useful without new matches."
    print("Evidence links from new stmts:", len(new_only_ev_links))
    print("Support links from new stmts:", len(new_only_mk_links))


    # Now merge in the "new" statements (the 20)
    updated_unique_stmts, updated_evidence_links, updated_match_key_links = \
        pas.merge_statements(original_unique_stmts, orignal_evidence_links,
                             original_match_key_links, new_stmts, optimize=True)

    # Check that we got all the same statements (trivial)
    assert len(updated_unique_stmts) == len(full_unique_stmts), \
        ("Got %d unique stmts from update, but %d from pre-assembly of all "
         "stmts." % (len(updated_unique_stmts), len(full_unique_stmts)))

    # Check that the evidence matches up (easy)
    fevl_set = pas.flatten_evidence_dict(full_evidence_links)
    uevl_set = pas.flatten_evidence_dict(updated_evidence_links)
    missed_evidence_links = fevl_set - uevl_set
    extra_evidence_links = uevl_set - fevl_set
    assert not len(missed_evidence_links) and not len(extra_evidence_links), \
        ("Some evidence links missed: %s, and/or some evidence links added: %s"
         % (_str_large_set(missed_evidence_links, 5),
            _str_large_set(extra_evidence_links, 5)))

    # Check that we have the same match key links (less easy)
    missed_match_keys = full_match_key_links - updated_match_key_links
    extra_match_keys = updated_match_key_links - full_match_key_links
    assert not len(extra_match_keys) and not len(missed_match_keys), \
        ("Some match key links missed: %s and/or some match key links added: %s"
         % (missed_match_keys, extra_match_keys))
