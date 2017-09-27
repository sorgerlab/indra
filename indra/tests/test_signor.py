from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.databases import hgnc_client
from indra.sources.signor import SignorProcessor, SignorRow, \
                                 _parse_residue_positions


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
    sp = SignorProcessor()
    assert isinstance(sp._data, list)
    assert isinstance(sp._data[0], SignorRow)
    assert isinstance(sp.statements, list)
    assert isinstance(sp.statements[0], Statement)


def test_get_agent():
    # Protein/gene
    test_ag = Agent('RELA', db_refs={'HGNC': hgnc_client.get_hgnc_id('RELA'),
                                     'UP': 'Q04206'})
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
    assert ev.source_api == 'signor'
    assert ev.source_id == 'SIGNOR-241929'
    assert ev.text == "Together, these results indicate that the Met gene " \
                      "is a direct target of NFkappaB and that Met " \
                      "participates in NFkappaB-mediated cell survival."


def test_process_row():
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)
    assert stmts[0].agent_list()[0].activity.activity_type == 'transcription'


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
    stmts, no_mech = SignorProcessor._process_row(test_row_chem_inh)
    assert not no_mech
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
    stmts, no_mech = SignorProcessor._process_row(test_row_chem_act)
    assert not no_mech
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
    stmts, no_mech = SignorProcessor._process_row(test_row_stab)
    assert not no_mech
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
    stmts, no_mech = SignorProcessor._process_row(test_row_destab)
    assert not no_mech
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
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
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
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].residue == 'T'
    assert stmts[1].position == '387'
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'T'
    assert mc.position == '387'
    assert af.is_active == True


def test_process_row_phos_down():
    test_row = SignorRow(ENTITYA='PRKCD', TYPEA='protein', IDA='Q05655',
            DATABASEA='UNIPROT', ENTITYB='PTPN22', TYPEB='protein',
            IDB='Q9Y2R2', DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='phosphorylation', RESIDUE='Ser35',
            SEQUENCE='FLKLKRQsTKYKADK', TAX_ID='9606', CELL_DATA='BTO:0000782',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='18056643', DIRECT='YES', NOTES='', ANNOTATOR='llicata',
            SENTENCE='', SIGNOR_ID='SIGNOR-159591')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Inhibition)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].residue == 'S'
    assert stmts[1].position == '35'
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'S'
    assert mc.position == '35'
    assert af.is_active == False


def test_process_row_phos_nores_up():
    test_row = SignorRow(ENTITYA='STK11', TYPEA='protein', IDA='Q15831',
            DATABASEA='UNIPROT', ENTITYB='AMPK', TYPEB='complex',
            IDB='SIGNOR-C15', DATABASEB='SIGNOR',
            EFFECT='up-regulates activity', MECHANISM='phosphorylation',
            RESIDUE='', SEQUENCE='', TAX_ID='-1', CELL_DATA='', TISSUE_DATA='', 
            MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='',
            MODASEQ='', MODIFICATIONB='', MODBSEQ='', PMID='14976552',
            DIRECT='YES', NOTES='', ANNOTATOR='lperfetto',
            SENTENCE='', SIGNOR_ID='SIGNOR-242602')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].residue is None
    assert stmts[1].position is None
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue is None
    assert mc.position is None
    assert af.is_active == True


def test_process_row_phos_nores_down():
    test_row = SignorRow(ENTITYA='CSNK1D', TYPEA='protein', IDA='P48730',
            DATABASEA='UNIPROT', ENTITYB='LEF1', TYPEB='protein', IDB='Q9UJU2',
            DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='phosphorylation', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='15747065', DIRECT='YES', NOTES='',
            ANNOTATOR='gcesareni', SENTENCE='', SIGNOR_ID='SIGNOR-134494')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Inhibition)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].residue is None
    assert stmts[1].position is None
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue is None
    assert mc.position is None
    assert af.is_active == False


def test_process_row_phos_multi_res():
    test_row = SignorRow(ENTITYA='RAF1', TYPEA='protein', IDA='P04049',
            DATABASEA='UNIPROT', ENTITYB='MAP2K2', TYPEB='protein',
            IDB='P36507', DATABASEB='UNIPROT', EFFECT='up-regulates',
            MECHANISM='phosphorylation', RESIDUE='Ser218;Ser222',
            SEQUENCE='VSGQLIDsMANSFVG;LIDSMANsFVGTRSY', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='8157000', DIRECT='YES',
            NOTES='', ANNOTATOR='gcesareni', SENTENCE='',
            SIGNOR_ID='SIGNOR-36553')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 4
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Phosphorylation)
    assert stmts[1].position == '218'
    assert isinstance(stmts[2], Phosphorylation)
    assert stmts[2].position == '222'
    af = stmts[3]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 2
    mc0 = af.agent.mods[0]
    assert mc0.mod_type == 'phosphorylation'
    assert mc0.residue == 'S'
    assert mc0.position == '218'
    mc1 = af.agent.mods[1]
    assert mc1.mod_type == 'phosphorylation'
    assert mc1.residue == 'S'
    assert mc1.position == '222'
    assert af.is_active == True


def test_process_row_complex_up():
    test_row = SignorRow(ENTITYA='NONO', TYPEA='protein', IDA='Q15233',
            DATABASEA='UNIPROT', ENTITYB='TOP1', TYPEB='protein', IDB='P11387',
            DATABASEB='UNIPROT', EFFECT='up-regulates', MECHANISM='binding',
            RESIDUE='', SEQUENCE='', TAX_ID='9606', CELL_DATA='BTO:0000017',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='9756848', DIRECT='YES', NOTES='', ANNOTATOR='miannu',
            SENTENCE='', SIGNOR_ID='SIGNOR-60557')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Complex)
    cplx_agent_a = stmts[1].agent_list()[0]
    cplx_agent_b = stmts[1].agent_list()[1]
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    # Won't fully match because of bound condition, so we check name
    assert af.agent.name == cplx_agent_b.name
    assert len(af.agent.bound_conditions) == 1
    bc = af.agent.bound_conditions[0]
    assert bc.agent.matches(cplx_agent_a)
    assert bc.is_bound
    assert af.activity == 'activity'
    assert af.is_active


def test_process_row_complex_down():
    test_row = SignorRow(ENTITYA='XIAP', TYPEA='protein', IDA='P98170',
            DATABASEA='UNIPROT', ENTITYB='CASP3', TYPEB='protein',
            IDB='P42574', DATABASEB='UNIPROT', EFFECT='down-regulates activity',
            MECHANISM='binding', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='10548111', DIRECT='YES', NOTES='',
            ANNOTATOR='amattioni', SENTENCE='', SIGNOR_ID='SIGNOR-71954')
    stmts, no_mech = SignorProcessor._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Inhibition)
    assert isinstance(stmts[1], Complex)
    cplx_agent_a = stmts[1].agent_list()[0]
    cplx_agent_b = stmts[1].agent_list()[1]
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    # Won't fully match because of bound condition, so we check name
    assert af.agent.name == cplx_agent_b.name
    assert len(af.agent.bound_conditions) == 1
    bc = af.agent.bound_conditions[0]
    assert bc.agent.matches(cplx_agent_a)
    assert bc.is_bound
    assert af.activity == 'activity'
    assert not af.is_active


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

