from __future__ import absolute_import, print_function, unicode_literals

import pytest
from indra.databases.lincs_client import get_drug_target_data, LincsClient


lc = LincsClient()


@pytest.mark.webservice
@unittest.skip('LINCS web service very unreliable.')
def test_get_drug_target_data():
    data_list = get_drug_target_data()
    assert len(data_list) > 100, len(data_list)


def test_get_protein_refs():
    prot_refs = lc.get_protein_refs('200020')
    assert prot_refs.get('UP') == 'P00519'
    assert prot_refs.get('EGID') == '25'
    assert prot_refs.get('HMS-LINCS') == '200020'


def test_get_sm_name():
    sm_name = lc.get_small_molecule_name('10001')
    assert sm_name == 'Seliciclib', sm_name


def test_get_sm_refs():
    sm_refs = lc.get_small_molecule_refs('10001')
    assert sm_refs.get('HMS-LINCS') == '10001', sm_refs
    assert sm_refs.get('PUBCHEM') == '160355', sm_refs
    assert sm_refs.get('CHEMBL') == 'CHEMBL14762', sm_refs

    sm_refs = lc.get_small_molecule_refs('10001-101')
    assert sm_refs.get('HMS-LINCS') == '10001-101', sm_refs
    assert sm_refs.get('PUBCHEM') == '160355', sm_refs
    assert sm_refs.get('CHEMBL') == 'CHEMBL14762', sm_refs
