from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.literature import pmc_client
from nose.tools import raises
from indra.util import unicode_strs

example_ids = {'pmid': '25361007',
               'pmcid': 'PMC4322985',
               'doi': '10.18632/oncotarget.2555'}

def test_id_lookup_pmid_no_prefix_no_idtype():
    ids = pmc_client.id_lookup('25361007')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_pmid_with_prefix_no_idtype():
    ids = pmc_client.id_lookup('PMID25361007')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_pmcid_no_idtype():
    ids = pmc_client.id_lookup('PMC4322985')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_pmcid_idtype():
    ids = pmc_client.id_lookup('PMC4322985', idtype='pmcid')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_pmcid_no_prefix_idtype():
    ids = pmc_client.id_lookup('4322985', idtype='pmcid')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_doi_no_prefix_no_idtype():
    ids = pmc_client.id_lookup('10.18632/oncotarget.2555')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

def test_id_lookup_doi_prefix_no_idtype():
    ids = pmc_client.id_lookup('DOI10.18632/oncotarget.2555')
    assert(ids['doi'] == example_ids['doi'])
    assert(ids['pmid'] == example_ids['pmid'])
    assert(ids['pmcid'] == example_ids['pmcid'])
    assert unicode_strs(ids)

@raises(ValueError)
def test_invalid_idtype():
    ids = pmc_client.id_lookup('DOI10.18632/oncotarget.2555', idtype='foo')

def test_get_xml():
    pmc_id = '4322985'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is not None)
    assert unicode_strs((pmc_id, xml_str))

def test_get_xml_PMC():
    pmc_id = 'PMC4322985'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is not None)
    assert unicode_strs((pmc_id, xml_str))

def test_get_xml_invalid():
    pmc_id = '9999999'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is None)

