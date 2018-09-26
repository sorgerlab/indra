from __future__ import absolute_import, print_function, unicode_literals

from indra.databases.lincs_client import get_drug_target_data, \
    get_small_molecule_data, get_protein_data, LincsClient


def test_get_drug_target_data():
    data_list = get_drug_target_data()
    assert len(data_list) > 100, len(data_list)


def test_get_small_molecule_data():
    sm_data = get_small_molecule_data()
    assert len(sm_data) > 100, len(sm_data)


def test_get_protein_data():
    prot_data = get_protein_data()
    assert len(prot_data) > 100, len(prot_data)


def test_lincs_client():
    lc = LincsClient()
    prot_refs = lc.get_protein_ref('25', id_type='entrez')
    assert prot_refs
    assert '200020' in prot_refs.keys()
    prot_ref_2 = lc.get_protein_ref('200020')  # id_type='hms-lincs'
    assert prot_refs['200020'] == prot_ref_2
    sm_name = lc.get_small_molecule_name('10001', id_type='short-hms-lincs')
    assert sm_name
    assert len(sm_name) == 1
    assert list(sm_name.values())[0] == 'Seliciclib'
    sm_refs = lc.get_small_molecule_ref('10001', id_type='short-hms-lincs')
    assert sm_refs
    assert set(sm_name.keys()) == set(sm_refs.keys())
    sm_ref = lc.get_small_molecule_ref('10001-101')  # id_type='hms-lincs'
    assert sm_ref == list(sm_refs.values())[0]
    sm_refs_2 = lc.get_small_molecule_ref('14762', id_type='chembl')
    assert sm_refs_2 == sm_refs
