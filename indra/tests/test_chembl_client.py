from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import Agent
from indra.databases import chembl_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr
import unittest

vem = Agent('VEMURAFENIB', db_refs={'CHEBI': '63637', 'TEXT': 'VEMURAFENIB'})
az628 = Agent('AZ628', db_refs={'CHEBI': '91354'})
braf = Agent('BRAF', db_refs={'HGNC': '1097', 'NCIT': 'C51194',
                              'TEXT': 'BRAF', 'UP': 'P15056'})
vem_chembl_id = 'CHEMBL1229517'
braf_chembl_id = 'CHEMBL5145'
query_dict_vem_activity = {'query': 'activity',
                           'params': {'molecule_chembl_id': vem_chembl_id,
                                      'limit': 10000}}
query_dict_BRAF_target = {'query': 'target',
                          'params': {'target_chembl_id': braf_chembl_id,
                                     'limit': 1}}


@attr('webservice')
def test_get_inhibitions():
    stmt = chembl_client.get_inhibition(vem, braf)
    assert(stmt is not None)
    assert(unicode_strs(stmt))
    assert(len(stmt.evidence) > 5)
    for ev in stmt.evidence:
        assert(ev.pmid)
        assert(ev.annotations)
        assert(ev.source_api == 'chembl')
        assert(ev.source_id)


@attr('webservice')
def test_activity_query():
    res = chembl_client.send_query(query_dict_vem_activity)
    assert(res['page_meta']['total_count'] == len(res['activities']))
    assay_types = list(set([x['standard_type'] for x in res['activities']]))
    expected_types = ['IC50', 'EC50', 'INH', 'Potency']
    for e_t in expected_types:
        assert(e_t in assay_types)


@attr('webservice')
def test_target_query():
    target = chembl_client.query_target(braf_chembl_id)
    assert(target['target_type'] == 'SINGLE PROTEIN')


@attr('webservice', 'slow')
def test_get_drug_inhibition_stmts_vem():
    stmts = chembl_client.get_drug_inhibition_stmts(vem)
    assert(len(stmts) > 0)
    for st in stmts:
        assert(unicode_strs(st))
        assert(len(st.evidence) >= 1)
        for ev in st.evidence:
            assert(ev.pmid)
            assert(ev.annotations)
            assert(ev.source_api == 'chembl')
            assert(ev.source_id)


@attr('webservice', 'slow')
def test_get_drug_inhibition_stmts_az628():
    stmts = chembl_client.get_drug_inhibition_stmts(az628)
    assert(len(stmts) > 0)
    for st in stmts:
        assert(unicode_strs(st))
        assert(len(st.evidence) >= 1)
        for ev in st.evidence:
            assert(ev.pmid)
            assert(ev.annotations)
            assert(ev.source_api == 'chembl')
            assert(ev.source_id)
