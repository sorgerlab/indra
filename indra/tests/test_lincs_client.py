from __future__ import absolute_import, print_function, unicode_literals

from indra.databases.lincs_client import get_drug_target_data,\
    get_small_molecule_data, get_protein_data


def test_get_drug_target_data():
    data_list = get_drug_target_data()
    assert len(data_list) > 100, len(data_list)


def test_get_small_molecule_data():
    sm_data = get_small_molecule_data()
    assert len(sm_data) > 100, len(sm_data)


def test_get_protein_data():
    prot_data = get_protein_data()
    assert len(prot_data) > 100, len(prot_data)
