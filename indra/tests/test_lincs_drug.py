from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.lincs.lincs_client import get_drug_target_data, \
    get_small_molecule_data, get_protein_data
from indra.sources.lincs.api import process_from_web


def test_get_drug_target_data():
    data_list = get_drug_target_data()
    assert len(data_list) > 100, len(data_list)


def test_get_small_molecule_data():
    sm_data = get_small_molecule_data()
    assert len(sm_data) > 100, len(sm_data)


def test_get_protein_data():
    prot_data = get_protein_data()
    assert len(prot_data) > 100, len(prot_data)


def test_process_from_web():
    lincs_p = process_from_web()
    assert lincs_p is not None
    assert lincs_p.statements
    data_len = len(get_drug_target_data())
    num_stmts = len(lincs_p.statements)
    assert num_stmts == data_len, \
        ("Did not convert all statements: expected %d, got %d."
         % (data_len, num_stmts))
