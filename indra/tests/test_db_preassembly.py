from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
import pickle
import random

from indra.db import util as db_util
from indra.db import preassembly_script as pas
from indra.statements import stmts_from_json


STMT_PICKLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'test_stmts_tuples.pkl')


STMTS = None
def _get_stmts():
    global STMTS
    if STMTS is None:
        with open(STMT_PICKLE_FILE, 'rb') as f:
            stmt_tuples = pickle.load(f)
        col_names = stmt_tuples.pop(0)
        stmt_jsons = [json.loads(tpl[-1].decode('utf8')) for tpl in stmt_tuples]
        STMTS = stmts_from_json(stmt_jsons)
    return STMTS


def test_preassembly_without_database():
    stmts = _get_stmts()
    unique_stmt_dict, evidence_links, match_key_links = \
        pas.process_statements(stmts)
    assert len(unique_stmt_dict)
    assert len(evidence_links) == len(stmts)
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
    _, _, new_only_match_key_links = pas.process_statements(new_stmts)
    assert len(new_only_match_key_links), "Test not useful without new matches."

    # Now merge in the "new" statements (the 20)
    updated_unique_stmts, updated_evidence_links, updated_match_key_links = \
        pas.merge_statements(original_unique_stmts, orignal_evidence_links,
                             original_match_key_links, new_stmts, optimize=True)

    # Check the results
    assert len(updated_unique_stmts) == len(full_unique_stmts)
    missed_evidence_links = full_evidence_links - updated_evidence_links
    extra_evidence_links = updated_evidence_links - full_evidence_links
    assert not len(missed_evidence_links) and not len(extra_evidence_links), \
        ("Some evidence links missed: %s, and/or some evidence links added: %s"
         % (missed_evidence_links, extra_evidence_links))
    missed_match_keys = full_match_key_links - updated_match_key_links
    extra_match_keys = updated_match_key_links - full_match_key_links
    assert not len(extra_match_keys) and not len(missed_match_keys), \
        ("Some match key links missed: %s and/or some match key links added: %s"
         % (missed_match_keys, extra_match_keys))
