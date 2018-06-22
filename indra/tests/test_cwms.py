from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises
from nose.plugins.attrib import attr

from indra.statements import *
from indra.sources.cwms import process_rdf_file, process_text

# Path to the CWMS test/dummy data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this, 'cwms_tests_data')

example1_rdf = os.path.join(data_folder, 'example_2_sentence_1.rdf')
example1_txt = os.path.join(data_folder, 'example_2_sentence_1.txt')

example2_rdf = os.path.join(data_folder, 'example_2_sentence_3.rdf')
example2_txt = os.path.join(data_folder, 'example_2_sentence_3.txt')

example3_rdf = os.path.join(data_folder, 'example_2_sentence_4.rdf')
example3_txt = os.path.join(data_folder, 'example_2_sentence_4.txt')


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
    subj = s0.subj
    assert subj.db_refs['TEXT'] == 'government', subj.db_refs['TEXT']
    assert subj.db_refs['CWMS'] == 'ONT::FEDERAL-ORGANIZATION',\
        subj.db_refs['CWMS']

    obj = s0.obj
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
    assert(isinstance(s0, Influence))
    subj = s0.subj
    assert(subj.db_refs['TEXT'] == 'Persistent insecurity and armed conflict')

    obj = s0.obj
    assert(obj.db_refs['TEXT'] == 'livelihood activities')

    ev = s0.evidence[0]
    assert(ev.text == text)
    assert(ev.source_api == 'cwms')


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
    assert subj.db_refs['TEXT'] == 'government', subj.db_refs['TEXT']
    assert subj.db_refs['CWMS'] == 'ONT::FEDERAL-ORGANIZATION',\
        subj.db_refs['CWMS']

    obj = s0.obj
    assert obj.db_refs['TEXT'] == 'agriculture', obj.db_refs['TEXT']
    assert obj.db_refs['CWMS'] == 'ONT::AGRICULTURE',\
        obj.db_refs['CWMS']

    ev = s0.evidence[0]
    assert ev.text == 'government influences agriculture.', ev.text
    assert ev.source_api == 'cwms', ev.source_api


def test_rdf_example1():
    # Example: These impacts on livestock and crops have resulted in
    # livelihoods being decimated.

    txt = load_text(example1_txt)
    cp = process_rdf_file(txt, example1_rdf)
    assert(len(cp.statements) == 1)

    statement0 = cp.statements[0]
    assert(statement0.subj.db_refs['TEXT'] ==
           'These impacts on livestock and crops')
    assert(statement0.obj.db_refs['TEXT'] ==
           'in livelihoods being decimated')


def test_rdf_example2():
    # Conflict and economic decline have led to violence and displacement.

    txt = load_text(example2_txt)
    cp = process_rdf_file(txt, example2_rdf)
    assert(len(cp.statements) == 1)

    statement0 = cp.statements[0]
    assert(statement0.subj.db_refs['TEXT'] ==
           'Conflict and economic decline')
    assert(statement0.obj.db_refs['TEXT'] ==
           'to violence and displacement')


def test_rdf_example3():
    # Violence has caused livestock to be looted, killed and disease-prone and
    # crops destroyed, and displacement has caused delayed planting.

    txt = load_text(example3_txt)
    cp = process_rdf_file(txt, example3_rdf)
    assert(len(cp.statements) == 1)

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
    text = 'Drought increases the decrease of crops by army worms'
    cp = process_text(text)
    assert cp is not None
    print(cp.statements)  # Check to make sure str/repr work.
    assert len(cp.statements) == 2, len(cp.statements)


@attr('slow', 'webservice')
def test_three_sentences():
    text = 'Floods cause displacement. Displacement reduces access to food. ' \
           'Rainfall causes floods.'
    cp = process_text(text)
    assert cp is not None
    print(cp.statements)
    assert len(cp.statements) == 3, len(cp.statements)
