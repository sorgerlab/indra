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
    stmts = SignorProcessor._process_row(test_row)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)
    assert stmts[0].agent_list()[0].activity.activity_type == 'transcription'


def test_process_row_binding():
    test_row_binding = SignorRow(ENTITYA='SERPINA1', TYPEA='protein',
        IDA='P01009', DATABASEA='UNIPROT', ENTITYB='LRP1', TYPEB='protein',
        IDB='Q07954', DATABASEB='UNIPROT', EFFECT='up-regulates',
        MECHANISM='binding', RESIDUE='', SEQUENCE='', TAX_ID='9606',
        CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
        MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
        PMID='8626456', DIRECT='YES', NOTES='', ANNOTATOR='gcesareni',
        SENTENCE='', SIGNOR_ID='SIGNOR-41180')
    stmts = SignorProcessor._process_row(test_row_binding)
    assert isinstance(stmts, list)
    assert len(stmts) == 2
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Complex)


def test_process_row_chem_inh():
    test_row_chem_inh = SignorRow(ENTITYA='722544-51-6', TYPEA='chemical',
        IDA='CID:16007391', DATABASEA='PUBCHEM', ENTITYB='AURKB',
        TYPEB='protein', IDB='Q96GD4', DATABASEB='UNIPROT',
        EFFECT='down-regulates', MECHANISM='chemical inhibition', RESIDUE='',
        SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='Other', DIRECT='YES',
        NOTES='Selleck', ANNOTATOR='gcesareni', SENTENCE='',
        SIGNOR_ID='SIGNOR-190245')
    stmts = SignorProcessor._process_row(test_row_chem_inh)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Inhibition)


def test_process_row_chem_act():
    test_row_chem_act = SignorRow(ENTITYA='Prostaglandin E2',
        TYPEA='smallmolecule', IDA='CID:5280360', DATABASEA='PUBCHEM',
        ENTITYB='GNG12', TYPEB='protein', IDB='Q9UBI6', DATABASEB='UNIPROT',
        EFFECT='up-regulates', MECHANISM='chemical activation', RESIDUE='',
        SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='16293724', DIRECT='YES', NOTES='',
        ANNOTATOR='gcesareni', SENTENCE='', SIGNOR_ID='SIGNOR-141820')
    stmts = SignorProcessor._process_row(test_row_chem_act)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Activation)


def test_process_row_stab():
    test_row_stab = SignorRow(ENTITYA='UCHL5', TYPEA='protein', IDA='Q9Y5K5',
            DATABASEA='UNIPROT', ENTITYB='TGFBR1', TYPEB='protein',
            IDB='P36897', DATABASEB='UNIPROT', EFFECT='up-regulates',
            MECHANISM='stabilization', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='17052192', DIRECT='YES', NOTES='',
            ANNOTATOR='gcesareni', SENTENCE='', SIGNOR_ID='SIGNOR-150135')
    stmts = SignorProcessor._process_row(test_row_stab)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)


def test_process_row_destab():
    test_row_destab = SignorRow(ENTITYA='INS', TYPEA='protein', IDA='P01308',
            DATABASEA='UNIPROT', ENTITYB='APOB', TYPEB='protein', IDB='P04114',
            DATABASEB='UNIPROT',
            EFFECT='down-regulates quantity by destabilization',
            MECHANISM='destabilization', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='BTO:0000575', TISSUE_DATA='',
            MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='',
            MODASEQ='', MODIFICATIONB='', MODBSEQ='', PMID='23721961',
            DIRECT='NO', NOTES='', ANNOTATOR='miannu',
            SENTENCE='', SIGNOR_ID='SIGNOR-252114')
    stmts = SignorProcessor._process_row(test_row_destab)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], DecreaseAmount)


def test_process_row_binding_complex():
    test_row = SignorRow(ENTITYA='ATG5', TYPEA='protein', IDA='Q9H1Y0',
            DATABASEA='UNIPROT', ENTITYB='ATG12/5/16L1', TYPEB='complex',
            IDB='SIGNOR-C109', DATABASEB='SIGNOR', EFFECT='form complex',
            MECHANISM='binding', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='BTO:0000007', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='18321988', DIRECT='YES', NOTES='',
            ANNOTATOR='lperfetto', SENTENCE='', SIGNOR_ID='SIGNOR-226693')
    stmts = SignorProcessor._process_row(test_row)
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Complex)
    assert len(stmts[0].agent_list()) == 2


def test_process_row_phos_up():
    test_row = SignorRow(ENTITYA='CHEK2', TYPEA='protein', IDA='O96017',
            DATABASEA='UNIPROT', ENTITYB='CHEK2', TYPEB='protein', IDB='O96017',
            DATABASEB='UNIPROT', EFFECT='up-regulates activity',
            MECHANISM='phosphorylation', RESIDUE='Thr387',
            SEQUENCE='LMRTLCGtPTYLAPE', TAX_ID='9606', CELL_DATA='BTO:0000007',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='11901158', DIRECT='YES', NOTES='', ANNOTATOR='gcesareni',
            SENTENCE='', SIGNOR_ID='SIGNOR-116131')
    stmts = SignorProcessor._process_row(test_row)
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].residue == 'T'
    assert stmts[1].position == '387'
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = ag.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'T'
    assert mc.position == '387'


def test_process_row_phos_down():
    # TODO
    pass


def test_process_row_complex_up():
    # TODO
    pass


def test_process_row_complex_down():
    # TODO
    pass


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
    #test_process_row_destab()
    #test_parse_csv()
    #test_get_agent()
    #test_get_agent_keyerror()
    #test_get_evidence()
    #test_process_row()
    #test_process_row_binding()
    #test_process_row_chem_inh()
    #test_parse_residue_positions()
    test_get_mechanism()
