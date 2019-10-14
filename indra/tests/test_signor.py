from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from os.path import join, dirname
from nose.tools import raises
from nose.plugins.attrib import attr

from indra.statements import *
from indra.databases import hgnc_client
from indra.sources.signor.processor import SignorProcessor, \
                                           _parse_residue_positions
from indra.sources.signor.api import _SignorRow_, process_from_file, \
                                     process_from_web

test_row = _SignorRow_(ENTITYA='RELA', TYPEA='protein', IDA='Q04206',
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
        SCORE='0.5',
        SIGNOR_ID='SIGNOR-241929')


test_data_file = join(dirname(__file__), 'signor_test_data.csv')
test_complexes_file = join(dirname(__file__), 'signor_test_complexes.csv')


def test_parse_csv_from_file():
    # Should work with both data file and complexes
    sp = process_from_file(test_data_file, test_complexes_file)
    assert isinstance(sp._data, list)
    assert isinstance(sp._data[0], _SignorRow_)
    assert isinstance(sp.statements, list)
    assert isinstance(sp.statements[0], Statement)
    # Test the complex map
    assert isinstance(sp.complex_map, dict)
    assert len(sp.complex_map) == 9
    assert 'SIGNOR-C1' in sp.complex_map
    assert isinstance(sp.complex_map['SIGNOR-C1'], list)
    assert sp.complex_map['SIGNOR-C1'] == ['P23511', 'P25208', 'Q13952']
    # Make sure we don't error if Complexes data is not provided
    sp = process_from_file(test_data_file)
    assert isinstance(sp.statements[0], Statement)
    assert sp.complex_map == {}


@attr('webservice', 'slow') 
def test_parse_csv_from_web():
    sp = process_from_web()
    assert isinstance(sp._data, list)
    assert isinstance(sp._data[0], _SignorRow_)
    assert isinstance(sp.statements, list)
    assert isinstance(sp.statements[0], Statement)
    # Test the complex map
    assert isinstance(sp.complex_map, dict)
    assert 'SIGNOR-C1' in sp.complex_map
    assert isinstance(sp.complex_map['SIGNOR-C1'], list)
    assert set(sp.complex_map['SIGNOR-C1']) == {'P23511', 'P25208', 'Q13952'}
    # Make sure we don't error if Complexes data is not provided
    for stmt in sp.statements:
        if isinstance(stmt, Complex):
            if len(stmt.members) < 2:
                assert False, 'Found a complex with less than 2 members: %s' %\
                    stmt


def test_get_agent():
    # Protein/gene
    # Create an empty Signor processor
    sp = SignorProcessor([])
    test_ag = Agent('RELA', db_refs={'HGNC': hgnc_client.get_hgnc_id('RELA'),
                                     'UP': 'Q04206'})
    sp_ag = sp._get_agent(test_row.ENTITYA, test_row.TYPEA,
                                       test_row.IDA, test_row.DATABASEA)
    assert test_ag.matches(sp_ag)
    # Chemical
    test_ag = Agent('AZD1480', db_refs={'PUBCHEM': '16659841'})
    sp_ag = sp._get_agent('AZD1480', 'chemical', 'CID: 16659841',
                                       'PUBCHEM')
    assert test_ag.matches(sp_ag)
    # Signor phenotype
    test_ag = Agent('Cell cycle progr.', db_refs={'SIGNOR': 'SIGNOR-PH42'})
    sp_ag = sp._get_agent('Cell cycle progr.', 'phenotype',
                                       'SIGNOR-PH42', 'SIGNOR')
    assert test_ag.matches(sp_ag)
    # Ungrounded -- couldn't find a real example in the dataset
    test_ag = Agent('Foobar', db_refs={})
    sp_ag = sp._get_agent('Foobar', 'pathway', None, None)
    assert test_ag.matches(sp_ag)
    sp_ag = sp._get_agent('Foobar', 'antibody', None, None)
    assert test_ag.matches(sp_ag)


@raises(KeyError)
def test_get_agent_keyerror():
    # Create an empty Signor processor
    sp = SignorProcessor([])
    sp_ag = sp._get_agent('foo', 'bar', None, None)


def test_get_evidence():
    ev = SignorProcessor._get_evidence(test_row)
    assert isinstance(ev, Evidence)
    assert ev.pmid == '19530226'
    assert ev.annotations == {
            'SEQUENCE': None,
            'MODULATOR_COMPLEX': None,
            'TARGET_COMPLEX': None,
            'MODIFICATIONA': None,
            'MODASEQ': None,
            'MODIFICATIONB': None,
            'MODBSEQ': None,
            'NOTES': None,
            'ANNOTATOR': 'gcesareni',
        }
    assert ev.context.species.db_refs['TAXONOMY'] == '10090'
    assert ev.context.cell_type.db_refs['BTO'] == '0002895'
    assert ev.epistemics['direct']
    assert ev.source_api == 'signor'
    assert ev.source_id == 'SIGNOR-241929'
    assert ev.text == "Together, these results indicate that the Met gene " \
                      "is a direct target of NFkappaB and that Met " \
                      "participates in NFkappaB-mediated cell survival."


def test_process_row():
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)
    assert stmts[0].agent_list()[0].activity.activity_type == 'transcription'


def test_process_row_chem_inh():
    test_row_chem_inh = _SignorRow_(ENTITYA='722544-51-6', TYPEA='chemical',
        IDA='CID:16007391', DATABASEA='PUBCHEM', ENTITYB='AURKB',
        TYPEB='protein', IDB='Q96GD4', DATABASEB='UNIPROT',
        EFFECT='down-regulates', MECHANISM='chemical inhibition', RESIDUE='',
        SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='Other', DIRECT='YES',
        NOTES='Selleck', ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
        SIGNOR_ID='SIGNOR-190245')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row_chem_inh)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Inhibition)


def test_process_row_chem_act():
    test_row_chem_act = _SignorRow_(ENTITYA='Prostaglandin E2',
        TYPEA='smallmolecule', IDA='CID:5280360', DATABASEA='PUBCHEM',
        ENTITYB='GNG12', TYPEB='protein', IDB='Q9UBI6', DATABASEB='UNIPROT',
        EFFECT='up-regulates', MECHANISM='chemical activation', RESIDUE='',
        SEQUENCE='', TAX_ID='9606', CELL_DATA='', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='16293724', DIRECT='YES', NOTES='',
        ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
        SIGNOR_ID='SIGNOR-141820')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row_chem_act)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Activation)


def test_process_row_stab():
    test_row_stab = _SignorRow_(ENTITYA='UCHL5', TYPEA='protein', IDA='Q9Y5K5',
            DATABASEA='UNIPROT', ENTITYB='TGFBR1', TYPEB='protein',
            IDB='P36897', DATABASEB='UNIPROT', EFFECT='up-regulates',
            MECHANISM='stabilization', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='17052192', DIRECT='YES', NOTES='',
            ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-150135')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row_stab)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], IncreaseAmount)


def test_process_row_destab():
    test_row_destab = _SignorRow_(ENTITYA='INS', TYPEA='protein', IDA='P01308',
            DATABASEA='UNIPROT', ENTITYB='APOB', TYPEB='protein', IDB='P04114',
            DATABASEB='UNIPROT',
            EFFECT='down-regulates quantity by destabilization',
            MECHANISM='destabilization', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='BTO:0000575', TISSUE_DATA='',
            MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='',
            MODASEQ='', MODIFICATIONB='', MODBSEQ='', PMID='23721961',
            DIRECT='NO', NOTES='', ANNOTATOR='miannu',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-252114')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row_destab)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], DecreaseAmount)


def test_process_row_binding_complex():
    test_row = _SignorRow_(ENTITYA='ATG5', TYPEA='protein', IDA='Q9H1Y0',
            DATABASEA='UNIPROT', ENTITYB='ATG12/5/16L1', TYPEB='complex',
            IDB='SIGNOR-C109', DATABASEB='SIGNOR', EFFECT='form complex',
            MECHANISM='binding', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='BTO:0000007', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='18321988', DIRECT='YES', NOTES='',
            ANNOTATOR='lperfetto', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-226693')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Complex)
    assert len(stmts[0].agent_list()) == 2


def test_process_row_phos_up():
    test_row = _SignorRow_(ENTITYA='CHEK2', TYPEA='protein', IDA='O96017',
            DATABASEA='UNIPROT', ENTITYB='CHEK2', TYPEB='protein', IDB='O96017',
            DATABASEB='UNIPROT', EFFECT='up-regulates activity',
            MECHANISM='phosphorylation', RESIDUE='Thr387',
            SEQUENCE='LMRTLCGtPTYLAPE', TAX_ID='9606', CELL_DATA='BTO:0000007',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='11901158', DIRECT='YES', NOTES='', ANNOTATOR='gcesareni',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-116131')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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
    test_row = _SignorRow_(ENTITYA='PRKCD', TYPEA='protein', IDA='Q05655',
            DATABASEA='UNIPROT', ENTITYB='PTPN22', TYPEB='protein',
            IDB='Q9Y2R2', DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='phosphorylation', RESIDUE='Ser35',
            SEQUENCE='FLKLKRQsTKYKADK', TAX_ID='9606', CELL_DATA='BTO:0000782',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='18056643', DIRECT='YES', NOTES='', ANNOTATOR='llicata',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-159591')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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
    test_row = _SignorRow_(ENTITYA='STK11', TYPEA='protein', IDA='Q15831',
            DATABASEA='UNIPROT', ENTITYB='AMPK', TYPEB='complex',
            IDB='SIGNOR-C15', DATABASEB='SIGNOR',
            EFFECT='up-regulates activity', MECHANISM='phosphorylation',
            RESIDUE='', SEQUENCE='', TAX_ID='-1', CELL_DATA='', TISSUE_DATA='', 
            MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='',
            MODASEQ='', MODIFICATIONB='', MODBSEQ='', PMID='14976552',
            DIRECT='YES', NOTES='', ANNOTATOR='lperfetto',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-242602')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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
    test_row = _SignorRow_(ENTITYA='CSNK1D', TYPEA='protein', IDA='P48730',
            DATABASEA='UNIPROT', ENTITYB='LEF1', TYPEB='protein', IDB='Q9UJU2',
            DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='phosphorylation', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='15747065', DIRECT='YES', NOTES='',
            ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-134494')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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


def test_process_row_dephos_up():
    test_row = _SignorRow_(ENTITYA='CHEK2', TYPEA='protein', IDA='O96017',
            DATABASEA='UNIPROT', ENTITYB='CHEK2', TYPEB='protein', IDB='O96017',
            DATABASEB='UNIPROT', EFFECT='up-regulates activity',
            MECHANISM='dephosphorylation', RESIDUE='Thr387',
            SEQUENCE='LMRTLCGtPTYLAPE', TAX_ID='9606', CELL_DATA='BTO:0000007',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='11901158', DIRECT='YES', NOTES='', ANNOTATOR='gcesareni',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-116131')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Dephosphorylation)
    assert stmts[1].residue == 'T'
    assert stmts[1].position == '387'
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'T'
    assert mc.position == '387'
    assert af.is_active == False


def test_process_row_dephos_down():
    test_row = _SignorRow_(ENTITYA='PRKCD', TYPEA='protein', IDA='Q05655',
            DATABASEA='UNIPROT', ENTITYB='PTPN22', TYPEB='protein',
            IDB='Q9Y2R2', DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='dephosphorylation', RESIDUE='Ser35',
            SEQUENCE='FLKLKRQsTKYKADK', TAX_ID='9606', CELL_DATA='BTO:0000782',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='18056643', DIRECT='YES', NOTES='', ANNOTATOR='llicata',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-159591')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Inhibition)
    assert isinstance(stmts[1], Dephosphorylation)
    assert stmts[1].residue == 'S'
    assert stmts[1].position == '35'
    af = stmts[2]
    assert isinstance(af, ActiveForm)
    assert len(af.agent.mods) == 1
    mc = af.agent.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'S'
    assert mc.position == '35'
    assert af.is_active == True


def test_mod_unknown_effect():
    test_row = _SignorRow_(ENTITYA='JAK2', TYPEA='protein', IDA='O60674',
            DATABASEA='UNIPROT', ENTITYB='JAK2', TYPEB='protein', IDB='O60674',
            DATABASEB='UNIPROT', EFFECT='unknown', MECHANISM='phosphorylation',
            RESIDUE='Tyr1007', SEQUENCE='VLPQDKEyYKVKEPG', TAX_ID='-1',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='9111318', DIRECT='YES', NOTES='', ANNOTATOR='',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-251358')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 1
    assert isinstance(stmts[0], Phosphorylation)
    assert stmts[0].residue == 'Y'
    assert stmts[0].position == '1007'


def test_process_row_dephos_nores_up():
    test_row = _SignorRow_(ENTITYA='STK11', TYPEA='protein', IDA='Q15831',
            DATABASEA='UNIPROT', ENTITYB='AMPK', TYPEB='complex',
            IDB='SIGNOR-C15', DATABASEB='SIGNOR',
            EFFECT='up-regulates activity', MECHANISM='dephosphorylation',
            RESIDUE='', SEQUENCE='', TAX_ID='-1', CELL_DATA='', TISSUE_DATA='', 
            MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='',
            MODASEQ='', MODIFICATIONB='', MODBSEQ='', PMID='14976552',
            DIRECT='YES', NOTES='', ANNOTATOR='lperfetto',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-242602')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Activation)
    assert isinstance(stmts[1], Dephosphorylation)
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


def test_process_row_dephos_nores_down():
    test_row = _SignorRow_(ENTITYA='CSNK1D', TYPEA='protein', IDA='P48730',
            DATABASEA='UNIPROT', ENTITYB='LEF1', TYPEB='protein', IDB='Q9UJU2',
            DATABASEB='UNIPROT', EFFECT='down-regulates',
            MECHANISM='dephosphorylation', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='15747065', DIRECT='YES', NOTES='',
            ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-134494')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
    assert not no_mech
    assert isinstance(stmts, list)
    assert len(stmts) == 3
    assert isinstance(stmts[0], Inhibition)
    assert isinstance(stmts[1], Dephosphorylation)
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


def test_process_row_phos_multi_res():
    test_row = _SignorRow_(ENTITYA='RAF1', TYPEA='protein', IDA='P04049',
            DATABASEA='UNIPROT', ENTITYB='MAP2K2', TYPEB='protein',
            IDB='P36507', DATABASEB='UNIPROT', EFFECT='up-regulates',
            MECHANISM='phosphorylation', RESIDUE='Ser218;Ser222',
            SEQUENCE='VSGQLIDsMANSFVG;LIDSMANsFVGTRSY', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='8157000', DIRECT='YES',
            NOTES='', ANNOTATOR='gcesareni', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-36553')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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
    test_row = _SignorRow_(ENTITYA='NONO', TYPEA='protein', IDA='Q15233',
            DATABASEA='UNIPROT', ENTITYB='TOP1', TYPEB='protein', IDB='P11387',
            DATABASEB='UNIPROT', EFFECT='up-regulates', MECHANISM='binding',
            RESIDUE='', SEQUENCE='', TAX_ID='9606', CELL_DATA='BTO:0000017',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='9756848', DIRECT='YES', NOTES='', ANNOTATOR='miannu',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-60557')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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
    test_row = _SignorRow_(ENTITYA='XIAP', TYPEA='protein', IDA='P98170',
            DATABASEA='UNIPROT', ENTITYB='CASP3', TYPEB='protein',
            IDB='P42574', DATABASEB='UNIPROT', EFFECT='down-regulates activity',
            MECHANISM='binding', RESIDUE='', SEQUENCE='', TAX_ID='9606',
            CELL_DATA='', TISSUE_DATA='', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='10548111', DIRECT='YES', NOTES='',
            ANNOTATOR='amattioni', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-71954')
    # Create an empty Signor processor
    sp = SignorProcessor([])
    stmts, no_mech = sp._process_row(test_row)
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


def test_signor_family_famplex_mapping():
    test_row = _SignorRow_(ENTITYA='TLRs', TYPEA='proteinfamily',
            IDA='SIGNOR-PF20', DATABASEA='SIGNOR',
            ENTITYB='Interferon Production', TYPEB='phenotype',
            IDB='SIGNOR-PH16', DATABASEB='SIGNOR', EFFECT='up-regulates',
            MECHANISM='', RESIDUE='', SEQUENCE='', TAX_ID='9606', CELL_DATA='',
            TISSUE_DATA='', MODULATOR_COMPLEX='', TARGET_COMPLEX='',
            MODIFICATIONA='', MODASEQ='', MODIFICATIONB='', MODBSEQ='',
            PMID='20404851', DIRECT='NO', NOTES='', ANNOTATOR='lperfetto',
            SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-216310')
    complex_map = {}
    sp = SignorProcessor([test_row], complex_map)
    statements = sp.statements
    assert len(statements) == 1
    s0 = statements[0]
    assert s0.subj.db_refs['FPLX'] == 'TLR'
    assert s0.subj.db_refs['SIGNOR'] == 'SIGNOR-PF20'
    assert s0.subj.name == 'TLR'


def test_signor_complexes():
    test_row = _SignorRow_(ENTITYA='NFY',
        TYPEA='complex', IDA='SIGNOR-C1', DATABASEA='SIGNOR', ENTITYB='ID1',
        TYPEB='protein', IDB='P41134', DATABASEB='UNIPROT',
        EFFECT='up-regulates quantity by expression',
        MECHANISM='transcriptional activation', RESIDUE='', SEQUENCE='',
        TAX_ID='9606', CELL_DATA='BTO:0000972', TISSUE_DATA='',
        MODULATOR_COMPLEX='', TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='',
        MODIFICATIONB='', MODBSEQ='', PMID='18025157', DIRECT='NO', NOTES='',
        ANNOTATOR='', SENTENCE='', SCORE='0.5', SIGNOR_ID='SIGNOR-255746')
    complex_map = {'SIGNOR-C1': ['P23511', 'P25208', 'Q13952']}
    sp = SignorProcessor([test_row], complex_map)
    assert isinstance(sp.statements, list)
    assert len(sp.statements) == 2

    s0 = sp.statements[0]
    assert isinstance(s0, IncreaseAmount)
    assert s0.subj.db_refs['UP'] == 'P23511'
    assert s0.subj.db_refs['HGNC'] == '7804'
    assert s0.subj.name == 'NFYA'
    assert s0.subj.bound_conditions[0].agent.db_refs['UP'] == 'P25208'
    assert s0.subj.bound_conditions[0].agent.name == 'NFYB'
    assert s0.subj.bound_conditions[0].agent.db_refs['HGNC'] == '7805'
    assert s0.subj.bound_conditions[1].agent.db_refs['UP'] == 'Q13952'
    assert s0.subj.bound_conditions[1].agent.name == 'NFYC'
    assert s0.subj.bound_conditions[1].agent.db_refs['HGNC'] == '7806'

    s1 = sp.statements[1]
    assert isinstance(s1, Complex)
    assert len(s1.evidence) == 1
    assert s1.evidence[0].source_api == 'signor'
    assert s1.evidence[0].source_id == 'SIGNOR-C1'
    assert s1.evidence[0].text
    correct_up_ids = set(['P23511', 'Q13952', 'P25208'])
    correct_hgnc_ids = set(['7804', '7805', '7806'])
    correct_names = set(['NFYA', 'NFYC', 'NFYB'])
    actual_up_ids = set([m.db_refs['UP'] for m in s1.members])
    actual_hgnc_ids = set([m.db_refs['HGNC'] for m in s1.members])
    actual_names = set([m.name for m in s1.members])
    assert actual_up_ids == correct_up_ids
    assert actual_hgnc_ids == correct_hgnc_ids
    assert actual_names == correct_names


def test_recursive_complexes():
    test_row = _SignorRow_(ENTITYA='PAX7/MLL2 complex', TYPEA='complex',
            IDA='SIGNOR-C91', DATABASEA='SIGNOR', ENTITYB='MYF5',
            TYPEB='protein', IDB='P13349', DATABASEB='UNIPROT',
            EFFECT='up-regulates quantity by expression',
            MECHANISM='transcriptional regulation', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='BTO:0002314',
            TISSUE_DATA='BTO:0000887;BTO:0001103', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='22863532', DIRECT='NO', NOTES='',
            ANNOTATOR='miannu', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-198641')
    complex_map = {
            'SIGNOR-C87': ['P61964', 'Q9UBL3', 'Q9C005', 'Q15291'],
            'SIGNOR-C88': ['SIGNOR-C87', 'O14686'],
            'SIGNOR-C91': ['SIGNOR-C88', 'P23759']}
    sp = SignorProcessor([test_row], complex_map)
    assert isinstance(sp.statements, list)
    assert len(sp.statements) == 4

    s0 = sp.statements[0]
    assert isinstance(s0, IncreaseAmount)
    bc = s0.subj.bound_conditions
    assert bc[0].agent.db_refs['UP'] == 'O14686'
    assert bc[1].agent.db_refs['UP'] == 'P61964'
    assert bc[2].agent.db_refs['UP'] == 'Q9UBL3'
    assert bc[3].agent.db_refs['UP'] == 'Q9C005'
    assert bc[4].agent.db_refs['UP'] == 'Q15291'

    assert isinstance(sp.statements[1], Complex)
    correct_complex_ups_1 = {'P61964', 'Q9UBL3', 'Q9C005', 'Q15291'}
    actual_complex_ups_1 = {m.db_refs['UP'] for m in sp.statements[1].members}
    assert correct_complex_ups_1 == actual_complex_ups_1

    assert isinstance(sp.statements[2], Complex)
    correct_complex_ups_2 = {'O14686', 'P61964', 'Q9UBL3', 'Q9C005', 'Q15291'}
    actual_complex_ups_2 = {m.db_refs['UP'] for m in sp.statements[2].members}
    assert correct_complex_ups_2 == actual_complex_ups_2

    assert isinstance(sp.statements[3], Complex)
    correct_complex_ups_3 = {'P23759', 'O14686', 'P61964', 'Q9UBL3', 'Q9C005',
                             'Q15291'}
    actual_complex_ups_3 = {m.db_refs['UP'] for m in sp.statements[3].members}
    assert correct_complex_ups_3 == actual_complex_ups_3


def test_complexes_with_families():
    test_row = _SignorRow_(ENTITYA='MAPK1', TYPEA='protein',
            IDA='P27361', DATABASEA='UNIPROT', ENTITYB='CDO/JLP/P38',
            TYPEB='complex', IDB='SIGNOR-C22', DATABASEB='SIGNOR',
            EFFECT='up-regulates quantity by expression',
            MECHANISM='transcriptional regulation', RESIDUE='', SEQUENCE='',
            TAX_ID='9606', CELL_DATA='BTO:0002314',
            TISSUE_DATA='BTO:0000887;BTO:0001103', MODULATOR_COMPLEX='',
            TARGET_COMPLEX='', MODIFICATIONA='', MODASEQ='', MODIFICATIONB='',
            MODBSEQ='', PMID='22863532', DIRECT='NO', NOTES='',
            ANNOTATOR='miannu', SENTENCE='', SCORE='0.5',
            SIGNOR_ID='SIGNOR-198641')
    complex_map = {'SIGNOR-C22' : ['O60271', 'Q4KMG0', 'SIGNOR-PF14']}
    sp = SignorProcessor([test_row], complex_map)
    assert isinstance(sp.statements, list)
    assert len(sp.statements) == 2

    assert isinstance(sp.statements[0], IncreaseAmount)
    obj = sp.statements[0].obj
    assert obj.db_refs['UP'] == 'O60271'
    assert len(obj.bound_conditions) == 2
    assert obj.bound_conditions[0].agent.db_refs['UP'] == 'Q4KMG0'
    assert obj.bound_conditions[1].agent.db_refs['SIGNOR'] == 'SIGNOR-PF14'
    assert obj.bound_conditions[1].agent.db_refs['FPLX'] == 'ROBO'

    s1 = sp.statements[1]
    assert isinstance(s1, Complex)
    members = s1.members
    assert members[0].db_refs['UP'] == 'O60271'
    assert members[1].db_refs['UP'] == 'Q4KMG0'
    assert members[2].db_refs['SIGNOR'] == 'SIGNOR-PF14'
    assert members[2].db_refs['FPLX'] == 'ROBO'


def test_recursively_expand_complex_constituents():
    complex_map = {
            'SIGNOR-C87': ['P61964', 'Q9UBL3', 'Q9C005', 'Q15291'],
            'SIGNOR-C88': ['SIGNOR-C87', 'O14686'],
            'SIGNOR-C91': ['SIGNOR-C88', 'P23759']}
    sp = SignorProcessor([test_row], complex_map)
    constituents = sp._recursively_lookup_complex('SIGNOR-C91')
    assert constituents == ['P23759', 'O14686', 'P61964', 'Q9UBL3', 'Q9C005',
                            'Q15291']
