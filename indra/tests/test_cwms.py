from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import unittest
from os.path import join, dirname, abspath
from nose.plugins.attrib import attr

from indra.statements import *
from indra.sources.cwms import *

# Path to the CWMS test/dummy data folder
path_this = dirname(abspath(__file__))
data_folder = join(path_this, 'cwms_tests_data')

example1_rdf = join(data_folder, 'example_2_sentence_1.rdf')
example1_txt = join(data_folder, 'example_2_sentence_1.txt')

example2_rdf = join(data_folder, 'example_2_sentence_3.rdf')
example2_txt = join(data_folder, 'example_2_sentence_3.txt')

example3_rdf = join(data_folder, 'example_2_sentence_4.rdf')
example3_txt = join(data_folder, 'example_2_sentence_4.txt')

ekb_processing_test_file = join(data_folder, 'ekb_processing_test.ekb')


def load_text(fname):
    with open(fname, 'r') as f:
        return f.read()


@attr('slow', 'webservice')
def test_cwmsreader_cause():
    # Test extraction of causal relations from the cwms reader service
    text = 'government causes agriculture.'
    cp = process_text(text)
    statements = cp.statements
    assert len(statements) == 1, len(statements)

    s0 = statements[0]
    assert isinstance(s0, Influence), type(s0)
    subj = s0.subj.concept
    assert subj.db_refs['TEXT'] == 'government', subj.db_refs['TEXT']
    assert subj.db_refs['CWMS'] == 'ONT::FEDERAL-ORGANIZATION',\
        subj.db_refs['CWMS']

    obj = s0.obj.concept
    assert obj.db_refs['TEXT'] == 'agriculture', obj.db_refs['TEXT']
    assert obj.db_refs['CWMS'] == 'ONT::AGRICULTURE',\
        obj.db_refs['CWMS']

    ev = s0.evidence[0]
    assert ev.text == 'government causes agriculture.', ev.text
    assert ev.source_api == 'cwms', ev.source_api


@attr('slow', 'webservice')
def test_cwmsreader_inhibit():
    # Test extraction of inhibition relations from the cwms reader service
    text = 'Persistent insecurity and armed conflict have disrupted ' + \
        'livelihood activities.'
    cp = process_text(text)
    statements = cp.statements
    print(statements)
    assert len(statements) == 1, len(statements)

    s0 = statements[0]
    print('Statement:', s0)
    assert isinstance(s0, Influence)
    subj = s0.subj.concept
    assert subj.db_refs['TEXT'] == 'Persistent insecurity and armed conflict'

    obj = s0.obj.concept
    assert obj.db_refs['TEXT'] == 'livelihood activities'

    ev = s0.evidence[0]
    assert ev.text == text
    assert ev.source_api == 'cwms'


@attr('slow', 'webservice')
def test_cwmsreader_influence():
    # Test extraction of causal relations from the cwms reader service
    text = 'government influences agriculture.'
    cp = process_text(text)
    statements = cp.statements
    assert len(statements) == 1, len(statements)

    s0 = statements[0]
    assert isinstance(s0, Influence), type(s0)
    subj = s0.subj
    assert subj.concept.db_refs['TEXT'] == 'government', \
        subj.concept.db_refs['TEXT']
    assert subj.concept.db_refs['CWMS'] == 'ONT::FEDERAL-ORGANIZATION', \
        subj.concept.db_refs['CWMS']

    obj = s0.obj
    assert obj.concept.db_refs['TEXT'] == 'agriculture', \
        obj.concept.db_refs['TEXT']
    assert obj.concept.db_refs['CWMS'] == 'ONT::AGRICULTURE', \
        obj.concept.db_refs['CWMS']

    ev = s0.evidence[0]
    assert ev.text == 'government influences agriculture.', ev.text
    assert ev.source_api == 'cwms', ev.source_api


def test_rdf_example1():
    # Example: These impacts on livestock and crops have resulted in
    # livelihoods being decimated.

    txt = load_text(example1_txt)
    cp = process_rdf_file(txt, example1_rdf)
    assert len(cp.statements) == 1

    statement0 = cp.statements[0]
    assert statement0.subj.db_refs['TEXT'] == \
           'These impacts on livestock and crops'
    assert statement0.obj.db_refs['TEXT'] == \
           'in livelihoods being decimated'


def test_rdf_example2():
    # Conflict and economic decline have led to violence and displacement.

    txt = load_text(example2_txt)
    cp = process_rdf_file(txt, example2_rdf)
    assert len(cp.statements) == 1

    statement0 = cp.statements[0]
    assert statement0.subj.db_refs['TEXT'] == \
           'Conflict and economic decline'
    assert statement0.obj.db_refs['TEXT'] == \
           'to violence and displacement'


def test_rdf_example3():
    # Violence has caused livestock to be looted, killed and disease-prone and
    # crops destroyed, and displacement has caused delayed planting.

    txt = load_text(example3_txt)
    cp = process_rdf_file(txt, example3_rdf)
    assert len(cp.statements) == 1

    statement0 = cp.statements[0]
    assert statement0.subj.db_refs['TEXT'] == 'displacement',\
        statement0.subj.db_refs['TEXT']
    assert statement0.obj.db_refs['TEXT'] == 'delayed planting',\
        statement0.obj.db_refs['TEXT']


@attr('slow', 'webservice')
def test_cwms_agriculture_increases():
    text = 'Agriculture increases food security.'
    cp = process_text(text)
    assert cp
    assert len(cp.statements) == 1


@attr('slow', 'webservice')
def test_cwms_two_sentences():
    text = 'Floods decrease agriculture. Agriculture increases food security.'
    cp = process_text(text)
    assert cp is not None
    assert len(cp.statements) == 2


@attr('slow', 'webservice')
def test_second_order_statements():
    # NOTE: the second order statement feature is being developed elsewhere,
    # however this test should still pass as is.
    text = 'Drought increases the decrease of crops by army worms'
    cp = process_text(text)
    assert cp is not None
    print(cp.statements)  # Check to make sure str/repr work.
    assert len(cp.statements) == 2, len(cp.statements)


@attr('slow', 'webservice')
def test_three_sentences():
    # These sentences were used in the June 2018 WM East and West coast
    # hackathons for creating a simple test model constructed from all the
    # readers, and utilizing other components.
    text = 'Floods cause displacement. Displacement reduces access to food. ' \
           'Rainfall causes floods.'
    cp = process_text(text)
    assert cp is not None
    print(cp.statements)
    assert len(cp.statements) == 3, len(cp.statements)
    assert all(isinstance(st, Influence) for st in cp.statements), cp.statements


@attr('slow', 'webservice')
def test_context_influence_obj():
    text = 'Hunger causes displacement in 2018 in South Sudan.'
    cp = process_text(text)
    stmt = cp.statements[0]
    assert isinstance(stmt.obj, Migration)
    cont = stmt.obj.context
    assert cont is not None
    assert cont.time and cont.locations


@attr('slow', 'webservice')
def test_context_influence_subj():
    text = 'Hunger in 2018 in South Sudan causes displacement.'
    cp = process_text(text)
    stmt = cp.statements[0]
    cont = stmt.subj.context
    assert cont is not None
    assert cont.time and cont.geo_location, cont


@attr('slow', 'webservice')
def test_context_influence_subj_obj():
    text = 'Hunger in 2018 causes displacement in South Sudan.'
    cp = process_text(text)
    stmt = cp.statements[0]
    assert stmt.subj.context.time and stmt.obj.context.locations


def test_ekb_process():
    cp = process_ekb_file(ekb_processing_test_file)
    assert len(cp.statements) == 1


def test_process_increase_event_ekb():
    fname = join(data_folder, 'cwms_increase.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1
    stmt = cp.statements[0]
    assert isinstance(stmt, Event)
    assert stmt.delta.polarity == 1, stmt.delta
    assert stmt.concept.name == 'food insecurity', stmt.concept.name
    assert stmt.context, stmt.context
    assert len(stmt.evidence) == 1
    ev = stmt.evidence[0]
    assert ev.source_api == 'cwms'
    assert ev.context is None
    assert ev.text is not None


def test_process_cause_decrease_event_ekb():
    fname = join(data_folder, 'cause_decrease_event.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1, cp.statements
    stmt = cp.statements[0]
    assert isinstance(stmt, Influence), stmt
    assert stmt.obj.delta.polarity == -1, stmt.obj.delta


def test_process_cause_increase_event_ekb():
    fname = join(data_folder, 'cause_increase_event.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1, cp.statements
    stmt = cp.statements[0]
    assert isinstance(stmt, Influence), stmt
    assert stmt.obj.delta.polarity == 1, stmt.obj.delta


def test_process_correlation():
    fname = join(data_folder, 'association.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1, cp.statements
    stmt = cp.statements[0]
    assert isinstance(stmt, Association), stmt
    assert stmt.members[0].concept.db_refs['CWMS'] == 'ONT::PRECIPITATION'
    assert stmt.members[1].concept.db_refs['CWMS'] == 'ONT::FLOODING'
    assert stmt.overall_polarity() is None


def test_process_migration1():
    fname = join(data_folder, 'migration_sentence1.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1
    stmt = cp.statements[0]
    assert isinstance(stmt, Migration)
    assert stmt.concept.name.startswith('In Sudan'), stmt.concept.name
    assert len(stmt.context.locations) == 1
    assert isinstance(stmt.context.locations[0]['location'], RefContext)
    assert stmt.context.locations[0]['location'].name == "Sudan"
    assert stmt.context.locations[0]['role'] == "destination"
    assert isinstance(stmt.context.time, TimeContext)
    assert stmt.context.time.text == "the month of April"
    assert isinstance(stmt.delta, QuantitativeState)
    assert stmt.delta.value == 23000
    assert stmt.delta.unit == "absolute"
    assert stmt.delta.modifier == "MORE"


def test_process_migration2():
    fname = join(data_folder, 'migration_sentence2.ekb')
    cp = process_ekb_file(fname)
    assert len(cp.statements) == 1
    stmt = cp.statements[0]
    assert isinstance(stmt, Migration)
    assert stmt.concept.name.startswith('Since the'), stmt.concept.name
    assert len(stmt.context.locations) == 2
    assert isinstance(stmt.context.locations[0]['location'], RefContext)
    assert stmt.context.locations[0]['location'].name == "Ethiopia"
    assert stmt.context.locations[0]['role'] == "destination"
    assert isinstance(stmt.context.locations[1]['location'], RefContext)
    assert stmt.context.locations[1]['location'].name == "South Sudan"
    assert stmt.context.locations[1]['role'] == "origin"
    assert isinstance(stmt.context.time, TimeContext)
    assert stmt.context.time.text == "the beginning of September 2016"
    assert isinstance(stmt.delta, QuantitativeState)
    assert stmt.delta.value == 40000
    assert stmt.delta.unit == "absolute"
    assert stmt.delta.modifier == "less_than"
