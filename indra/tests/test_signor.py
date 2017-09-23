from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.databases import hgnc_client
from indra.sources.signor import SignorProcessor, SignorRow

def _id(gene):
    return hgnc_client.get_hgnc_id(gene)

signor_test_path = join(dirname(__file__), '..', '..', 'data',
                        'all_data_23_09_17.csv')

test_row = SignorRow(ENTITYA='RELA', TYPEA='protein', IDA='Q04206',
        DATABASEA='UNIPROT', ENTITYB='MET', TYPEB='protein', IDB='P08581',
        DATABASEB='UNIPROT', EFFECT='up-regulates quantity',
        MECHANISM='transcriptional regulation', RESIDUE='', SEQUENCE='',
        TAX_ID='10090', CELL_DATA='BTO:0002895', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='19530226', DIRECT='YES', NOTES='',
        ANNOTATOR='gcesareni',
        SENTENCE="Together, these results indicate that the Met gene is a "\
                 "direct target of NFkappaB and that Met participates "\
                 "in NFkappaB-mediated cell survival.",
        SIGNOR_ID='SIGNOR-241929')


def test_parse_csv():
    sp = SignorProcessor(signor_test_path)
    assert isinstance(sp._data, list)
    assert isinstance(sp._data[0], SignorRow)


def test_get_agent():
    # Protein/gene
    test_ag = Agent('RELA', db_refs={'HGNC': _id('RELA'), 'UP': 'Q04206'})
    sp_ag = SignorProcessor._get_agent(test_row.ENTITYA, test_row.TYPEA,
                                       test_row.IDA, test_row.DATABASEA)
    assert test_ag.matches(sp_ag)
    # Chemical
    test_ag = Agent('AZD1480', db_refs={'PUBCHEM': 'CID:16659841'})
    sp_ag = SignorProcessor._get_agent('AZD1480', 'chemical', 'CID:16659841',
                                       'PUBCHEM')
    assert test_ag.matches(sp_ag)
    # Signor phenotype
    test_ag = Agent('Cell cycle progr.', db_refs={'SIGNOR': 'SIGNOR-PH42'})
    sp_ag = SignorProcessor._get_agent('Cell cycle progr.', 'phenotype',
                                       'SIGNOR-PH42', 'SIGNOR')
    assert test_ag.matches(sp_ag)
    # Ungrounded -- couldn't find a real example in the dataset
    test_ag = Agent('Foobar', db_refs={})
    sp_ag = SignorProcessor._get_agent('Foobar', 'pathway', None, None)
    assert test_ag.matches(sp_ag)
    sp_ag = SignorProcessor._get_agent('Foobar', 'antibody', None, None)
    assert test_ag.matches(sp_ag)


@raises(KeyError)
def test_get_agent_keyerror():
    sp_ag = SignorProcessor._get_agent('foo', 'bar', None, None)


def test_get_evidence():
    ev = SignorProcessor._get_evidence(test_row)
    assert isinstance(ev, Evidence)
    assert ev.pmid == '19530226'
    assert ev.annotations == {
            'SEQUENCE': None,
            'TAX_ID': '10090',
            'CELL_DATA': 'BTO:0002895',
            'TISSUE_DATA': None,
            'MODULATOR_COMPLEX': None,
            'TARGET_COMPLEX': None,
            'MODIFICATIONA': None,
            'MODASEQ': None,
            'MODIFICATIONB': None,
            'MODBSEQ': None,
            'NOTES': None,
            'ANNOTATOR': 'gcesareni',
        }
    assert ev.epistemics['direct']
    assert ev.source_api == 'SIGNOR'
    assert ev.source_id == 'SIGNOR-241929'
    assert ev.text == "Together, these results indicate that the Met gene " \
                      "is a direct target of NFkappaB and that Met " \
                      "participates in NFkappaB-mediated cell survival."


def test_process_row():
    stmt = SignorProcessor._process_row(test_row)
    assert isinstance(stmt, IncreaseAmount)


def test_get_mechanism():
    sp = SignorProcessor(signor_test_path)
    assert sp.statements
    globals().update(locals())

if __name__ == '__main__':
    test_parse_csv()
    test_get_agent()
    test_get_agent_keyerror()
    test_get_evidence()
    test_process_row()
    test_get_mechanism()
