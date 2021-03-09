from indra.sources.minerva.api import *
from indra.sources.minerva.processor import SifProcessor
from indra.statements import Activation, Inhibition


model_id_to_sif_strs = {
    790: ['sa44 POSITIVE csa5', 'sa18 NEGATIVE csa3'],  # TGFB pathway model
    799: ['sa18 POSITIVE sa15', 'csa2 POSITIVE sa9']  # apoptosis pathway model
}


def test_process_sif_strs():
    sp = process_sif_strs(model_id_to_sif_strs)
    assert sp
    assert isinstance(sp, SifProcessor)
    assert len(sp.statements) == 4
    # Correct statement types
    assert isinstance(sp.statements[0], Activation)
    assert isinstance(sp.statements[1], Inhibition)
    assert isinstance(sp.statements[2], Activation)
    assert isinstance(sp.statements[3], Activation)
    # Using the same code ('sa18'), get different agents depending on model
    assert sp.statements[1].subj.name == 'MAPK3'
    assert sp.statements[1].subj.get_grounding() == ('HGNC', '6877')
    assert sp.statements[2].subj.name == 'CASP9'
    assert sp.statements[2].subj.get_grounding() == ('HGNC', '1511')
    # If possible, get FamPlex family
    # "sa44 POSITIVE csa5" is "ACVR1 POSITIVE SMAD2/3_complex"
    assert sp.statements[0].obj.name == 'SMAD2_3'
    assert sp.statements[0].obj.get_grounding() == ('FPLX', 'SMAD2_3')
    # Otherwise create Agent with BoundConditions
    # "csa2 POSITIVE sa9" is "csa2 POSITIVE sa9"
    assert sp.statements[3].subj.name == 'FAS'
    assert sp.statements[3].subj.bound_conditions
    assert sp.statements[3].subj.bound_conditions[0].agent.name == 'FASLG'
    # Both main agent and BoundCondition agent have groundings
    assert sp.statements[3].subj.get_grounding() == ('HGNC', '11920')
    assert sp.statements[3].subj.bound_conditions[0].agent.get_grounding() == (
        'HGNC', '11936')


def test_process_file():
    sp = process_file('minerva_test1.sif', 790)
    assert sp
    assert isinstance(sp, SifProcessor)
    assert len(sp.statements) == 2


def test_process_files():
    # One file
    sp = process_files({790: 'minerva_test1.sif'})
    assert sp
    assert isinstance(sp, SifProcessor)
    assert len(sp.statements) == 2
    # Multiple files
    sp = process_files({790: 'minerva_test1.sif', 799: 'minerva_test2.sif'})
    assert sp
    assert isinstance(sp, SifProcessor)
    assert len(sp.statements) == 4


def test_process_from_web():
    sp = process_from_web(model_ids=[790])
    assert sp
    assert isinstance(sp, SifProcessor)
    assert len(sp.statements) > 20
