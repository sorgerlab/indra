from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.databases import hgnc_client
from indra.sources.signor import SignorProcessor, SignorRow, \
                                 _parse_residue_positions

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

test_row_binding = SignorRow(ENTITYA='SERPINA1', TYPEA='protein', IDA='P01009', 
        DATABASEA='UNIPROT', ENTITYB='LRP1', TYPEB='protein', IDB='Q07954',
        DATABASEB='UNIPROT', EFFECT='up-regulates', MECHANISM='binding',
        RESIDUE='', SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='8626456', DIRECT='YES', NOTES='',
        ANNOTATOR='gcesareni', SENTENCE='In vitro binding studies revealed '\
                'that antithrombin iii (atiii)thrombin, heparin cofactor ii '\
                '(hcii)thrombin, and ?1-antitrypsin (?1AT)trypsin bound to '\
                'purified lrp',
        SIGNOR_ID='SIGNOR-41180')

test_row_dup1 = SignorRow(ENTITYA='722544-51-6', TYPEA='chemical',
        IDA='CID:16007391', DATABASEA='PUBCHEM', ENTITYB='AURKB',
        TYPEB='protein', IDB='Q96GD4', DATABASEB='UNIPROT',
        EFFECT='down-regulates', MECHANISM='chemical inhibition', RESIDUE='',
        SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='Other', DIRECT='YES',
        NOTES='Selleck', ANNOTATOR='gcesareni', SENTENCE='',
        SIGNOR_ID='SIGNOR-190245')

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
    (effect_stmt, mech_stmts, af_stmt) = SignorProcessor._process_row(test_row)
    assert isinstance(effect_stmt, IncreaseAmount)
    assert isinstance(mech_stmts, list)
    assert len(mech_stmts) == 0


def test_process_row_binding():
    (effect_stmt, mech_stmts, af_stmt) = \
                    SignorProcessor._process_row(test_row_binding)
    assert isinstance(effect_stmt, Activation)
    assert isinstance(mech_stmts, list)
    assert len(mech_stmts) == 1
    assert isinstance(mech_stmts[0], Complex)


def test_process_row_dup1():
    (effect_stmt, mech_stmts, af_stmt) = \
                    SignorProcessor._process_row(test_row_dup1)
    assert isinstance(effect_stmt, Inhibition)
    assert isinstance(mech_stmts, list)
    assert len(mech_stmts) == 0


def test_parse_residue_positions():
    residues = _parse_residue_positions('TYR304')
    assert len(residues) == 1
    assert residues[0][0] == 'Y'
    assert residues[0][1] == '304'
    # Invalid residue
    residues = _parse_residue_positions('Foo')
    assert len(residues) == 1
    assert residues[0] == (None, None)
    # Residue but not position
    residues = _parse_residue_positions('gly')
    assert residues[0][0] == 'G'
    assert residues[0][1] == None
    # Position can't be converted to int
    residues = _parse_residue_positions('glyxxx')
    assert len(residues) == 1
    assert residues[0] == (None, None)
    # Multiple positions separated by semicolons
    residues = _parse_residue_positions('Tyr1185; Tyr1190')
    assert len(residues) == 2
    assert residues[0][0] == 'Y'
    assert residues[0][1] == '1185'
    assert residues[1][0] == 'Y'
    assert residues[1][1] == '1190'
    residues = _parse_residue_positions('Thr169;Tyr171')
    assert len(residues) == 2
    assert residues[0][0] == 'T'
    assert residues[0][1] == '169'
    assert residues[1][0] == 'Y'
    assert residues[1][1] == '171'


def test_get_mechanism():
    sp = SignorProcessor(signor_test_path)
    assert sp.statements
    globals().update(locals())

if __name__ == '__main__':
    test_process_row_dup1()
    """
    test_parse_csv()
    test_get_agent()
    test_get_agent_keyerror()
    test_get_evidence()
    test_process_row()
    test_process_row_binding()
    test_parse_residue_positions()
    test_get_mechanism()
    """
