from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname

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
        ANNOTATOR='gcesareni', SENTENCE="""Together, these results indicate that
        the Met gene is a direct target of NFkappaB and that Met participates
        in NFkappaB-mediated cell survival.""", SIGNOR_ID='SIGNOR-241929')


def test_parse_csv():
    sp = SignorProcessor(signor_test_path)
    assert isinstance(sp._data, list)
    assert isinstance(sp._data[0], SignorRow)
    globals().update(locals())


def test_get_agent():
    # Protein/gene
    test_ag = Agent('RELA', db_refs={'HGNC': _id('RELA'), 'UP': 'Q04206'})
    sp_ag = sp._get_agent(test_row.ENTITYA, test_row.TYPEA, test_row.IDA,
                          test_row.DATABASEA)
    assert test_ag.matches(sp_ag)
    # Chemical
    test_ag = Agent('AZD1480', db_refs={'PUBCHEM': 'CID:16659841'})
    sp_ag = sp._get_agent('AZD1480', 'chemical', 'CID:16659841', 'PUBCHEM')
    assert test_ag.matches(sp_ag)


if __name__ == '__main__':
    test_parse_csv()
    test_get_agent()
