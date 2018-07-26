from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import unittest
from os.path import join, dirname
from indra.statements import *
from indra.sources.hume.api import *


# Path to the HUME test files
path_this = os.path.dirname(os.path.abspath(__file__))
test_file_simple = os.path.join(path_this, 'bbn_test_simple.json-ld')
test_file_negatedCause = os.path.join(path_this,
                                      'bbn_test_negatedCause.json-ld')
test_file_negatedEffect = os.path.join(path_this,
                                       'bbn_test_negatedEffect.json-ld')


def test_simple_extraction():
    """Verify that processor extracts a simple causal assertion correctly from
    a JSON-LD file."""
    bp = process_json_file_old(test_file_simple)
    statements = bp.statements

    assert(len(statements) == 1)
    s0 = statements[0]

    assert(isinstance(s0, Influence))
    assert(s0.subj.name == 'cow')
    assert(s0.subj.db_refs['HUME'] == 'Bovine')
    assert(s0.obj.name == 'moo')
    assert(s0.obj.db_refs['HUME'] == 'MooSound')

    assert(len(s0.evidence) == 1)
    ev0 = s0.evidence[0]
    assert(ev0.source_api == 'hume')
    assert(ev0.text == 'Cow causes moo.')


def test_negated_cause():
    """We only want to extract causal relations between two positive events.
    The processor should give no statements for a negated cause."""
    bp = process_json_file_old(test_file_negatedCause)
    assert(len(bp.statements) == 0)


def test_negated_effect():
    """We only want to extract causal relations between two positive events.
    The processor should give no statements for a negated effect."""
    bp = process_json_file_old(test_file_negatedEffect)
    assert(len(bp.statements) == 0)


@unittest.skip('Need updated JSON-LD file')
def test_bbn_on_ben_paragraph():
    bp = process_jsonld_file(join(path_this,
                                  'hackathon_test_paragraph.json-ld'))
    assert bp is not None
    print(bp.statements)
    stmt_dict = {hash(s.get_hash()): s for s in bp.statements}
    assert len(stmt_dict) == 3, len(stmt_dict)
