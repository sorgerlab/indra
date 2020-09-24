import os
import unittest
from builtins import dict
from datetime import datetime

from indra.sources.hume.api import *
from indra.statements import *

# Path to the HUME test files
path_this = os.path.dirname(os.path.abspath(__file__))

test_file_new_simple = os.path.join(path_this, 'wm_m12.ben_sentence.json-ld')

standalone_events = os.path.join(
    path_this, 'wm_ben_event_sentences.v1.json-ld')

migration_events = os.path.join(
    path_this, 'wm_migration_numeric_one_sentence.082019.json-ld')


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
    assert stmt.obj.context is not None
    assert stmt.obj.context.time is not None
    time = stmt.obj.context.time
    assert time.text == '2018', time.text
    assert isinstance(time.start, datetime), type(time.start)
    assert isinstance(time.end, datetime), type(time.end)
    assert isinstance(time.duration, int), type(time.duration)
    assert stmt.obj.context.geo_location is not None
    loc = stmt.obj.context.geo_location
    assert loc.name == 'South Sudan', loc.name


def test_standalone_events():
    bp = process_jsonld_file(standalone_events)
    assert bp, "Processor is none."
    assert len(bp.statements) == 3, len(bp.statements)
    food_stmt = [
        st for st in bp.statements if st.concept.name == 'insecurity'][0]
    conflict_stmt = [
        st for st in bp.statements if st.concept.name == 'Conflict'][0]
    assert isinstance(food_stmt, Event)
    assert food_stmt.context.geo_location.name == 'South Sudan', \
        food_stmt.context.geo_location.name
    assert food_stmt.context.time.text == '2019', food_stmt.context.time.text
    assert food_stmt.delta.polarity == 1
    assert len(food_stmt.evidence) == 1, len(food_stmt.evidence)
    assert isinstance(conflict_stmt, Event)
    assert conflict_stmt.context.time.text == 'May 2017', \
        conflict_stmt.context.time.text
    assert len(conflict_stmt.evidence) == 1, len(conflict_stmt.evidence)


def test_migration_events():
    bp = process_jsonld_file(migration_events, extract_filter={'event'})
    assert bp, "Processor is none."
    assert len(bp.statements) == 1
    stmt = bp.statements[0]
    assert isinstance(stmt, Migration)
    assert len(stmt.context.locations) == 2
    location_dict = dict()
    assert isinstance(stmt.context.locations[0]['location'], RefContext)
    for location_ref in stmt.context.locations:
        location_dict[location_ref['role']] = location_ref['location'].name
    assert location_dict['origin'] == "South Sudan"
    assert location_dict['destination'] == (
        "Federal Democratic Republic of Ethiopia")
    assert isinstance(stmt.context.time, TimeContext)
    assert "2018" in stmt.context.time.text
    assert isinstance(stmt.delta, QuantitativeState)
    assert stmt.delta.value == 10000
    assert stmt.delta.unit == "Monthly"
    assert stmt.delta.modifier == "Min"

    # Test extraction filter
    bp = process_jsonld_file(migration_events, extract_filter={'influence'})
    assert len(bp.statements) == 0


def test_compositional_grounding():
    fname = os.path.join(path_this, 'hume.compositional.output.json-ld')
    bp = process_jsonld_file(fname)
    assert bp
    assert bp.statements
    for stmt in bp.statements:
        for concept in stmt.agent_list():
            refs = concept.db_refs
            assert 'WM' in refs
            wm = refs['WM']
            assert len(wm) == 1, refs
            assert wm[0] is not None
            wmg = wm[0]
            assert all(len(entry) == 2 for entry in wmg if entry is not None)
            assert all(entry[0].startswith('wm_compositional') for entry
                       in wmg if entry is not None)
