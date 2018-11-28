from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import unittest
from datetime import datetime
from indra.statements import *
from indra.sources.hume.api import *

# Path to the HUME test files
path_this = os.path.dirname(os.path.abspath(__file__))

test_file_new_simple = os.path.join(path_this, 'wm_m12.ben_sentence.json-ld')


@unittest.skip('Need updated JSON-LD file')
def test_bbn_on_ben_paragraph():
    bp = process_jsonld_file(os.path.join(path_this,
                                          'hackathon_test_paragraph.json-ld'))
    assert bp is not None
    print(bp.statements)
    stmt_dict = {s.get_hash(shallow=False): s for s in bp.statements}
    assert len(stmt_dict) == 3, len(stmt_dict)


def test_large_bbn_corpus():
    file_path = os.path.join(path_this,
                             'wm_m12.v8.full.v4.json-ld')
    if not os.path.exists(file_path):
        raise unittest.SkipTest("The test file is not available.")
    bp = process_jsonld_file(os.path.join(path_this,
                             'wm_m12.v8.full.v4.json-ld'))
    assert bp is not None
    assert len(bp.statements) > 1000
    print(len(bp.statements))


def test_for_context():
    bp = process_jsonld_file(test_file_new_simple)
    assert bp, "Processor is none."
    assert len(bp.statements) == 1, len(bp.statements)
    stmt = bp.statements[0]
    assert len(stmt.evidence) == 1, len(stmt.evidence)
    ev = stmt.evidence[0]
    assert ev.context is not None
    assert ev.context.time is not None
    time = ev.context.time
    assert time.text == '2018', ev.context.time.text
    assert isinstance(time.start, datetime), type(time.start)
    assert isinstance(time.end, datetime), type(time.end)
    assert isinstance(time.duration, int), type(time.duration)
    assert ev.context.geo_location is not None
    loc = ev.context.geo_location
    assert loc.name == 'South Sudan', loc.name
