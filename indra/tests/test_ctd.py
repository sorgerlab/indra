import os
from indra.statements import *
from indra.sources import ctd
from indra.sources.ctd.processor import CTDChemicalGeneProcessor

HERE = os.path.dirname(os.path.abspath(__file__))


def test_statement_type_mapping():
    st = CTDChemicalGeneProcessor.get_statement_types(
        'decreases^phosphorylation', 'X',
        'X decreases the phosphorylation of Y')
    assert set(st.values()) == {Dephosphorylation}, st

    st = CTDChemicalGeneProcessor.get_statement_types(
        'decreases^reaction|increases^phosphorylation', 'X',
        'X decreases the reaction [Z increases the phosphorylation of Y]')
    assert set(st.values()) == {Dephosphorylation}, st


def test_chemical_gene():
    fname = os.path.join(HERE, 'ctd_chem_gene_20522546.tsv')
    cp = ctd.process_tsv(fname, 'chemical_gene')
    assert len(cp.statements) == 4, cp.statements
