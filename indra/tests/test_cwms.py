from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.sources.cwms import process_rdf_file

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

def test_rdf_example1():
    # Example: These impacts on livestock and crops have resulted in
    # livelihoods being decimated.

    txt = load_text(example1_txt)
    cp = process_rdf_file(txt, example1_rdf)
    assert(len(cp.statements) == 1)

    statement0 = cp.statements[0]
    assert(statement0.subj.db_refs['TEXT'] == \
            'These impacts on livestock and crops')
    assert(statement0.obj.db_refs['TEXT'] == \
            'in livelihoods being decimated')

def test_rdf_example2():
    # Conflict and economic decline have led to violence and displacement.

    txt = load_text(example2_txt)
    cp = process_rdf_file(txt, example2_rdf)
    assert(len(cp.statements) == 1)

    statement0 = cp.statements[0]
    assert(statement0.subj.db_refs['TEXT'] == \
            'Conflict and economic decline')
    assert(statement0.obj.db_refs['TEXT'] == \
            'to violence and displacement')

def test_rdf_example3():
    # Violence has caused livestock to be looted, killed and disease-prone and
    # crops destroyed, and displacement has caused delayed planting.
    
    txt = load_text(example3_txt)
    cp = process_rdf_file(txt, example3_rdf)
    assert(len(cp.statements) == 1)

    statement0 = cp.statements[0]
    assert(statement0.subj.db_refs['TEXT'] == \
            'displacement')
    assert(statement0.obj.db_refs['TEXT'] == \
            'delayed planting')

