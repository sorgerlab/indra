import pybel
from pybel.dsl import *
import pybel.constants as pc
from pybel.examples import egf_graph, sialic_acid_graph
from indra.statements import *
from indra.sources import pybel as pb
from indra.databases import hgnc_client
from nose.tools import raises, ok_

mek_hgnc_id = hgnc_client.get_hgnc_id('MAP2K1')
mek_up_id = hgnc_client.get_uniprot_id(mek_hgnc_id)


def test_process_pybel():
    pbp = pb.process_pybel_graph(egf_graph)
    assert pbp.statements


def test_get_agent_hgnc():
    mek = protein(name='MAP2K1', namespace='HGNC')
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id

    # Now create an agent with an identifier
    mek = protein(name='Foo', namespace='HGNC', identifier='6840')
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id


def test_get_agent_up():
    mek = protein(namespace='UP', identifier='Q02750')
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id


def test_get_agent_up_no_id():
    mek = protein(name='MAP2K1', namespace='UP')
    agent = pb._get_agent(mek, {})
    assert agent is None


def test_get_agent_with_mods():
    mek = protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph')])
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert not mod.residue
    assert not mod.position

    mek = protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', code='Ser')])
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert mod.residue == 'S'
    assert not mod.position

    mek = protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', position=218)])
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert not mod.residue
    assert mod.position == '218'

    mek = protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', position=218, code='Ser')])
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert mod.residue == 'S'
    assert mod.position == '218'


def test_get_agent_with_muts():
    mek = protein(name='MAP2K1', namespace='HGNC',
                  variants=[hgvs('p.Val600Glu')])
    agent = pb._get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mutations) == 1
    mut = agent.mutations[0]
    assert mut.position == '600'
    assert mut.residue_from == 'V'
    assert mut.residue_to == 'E'


def test_get_agent_with_activity():
    mek = protein(name='MAP2K1', namespace='HGNC')
    agent = pb._get_agent(mek, activity('act'))
    assert isinstance(agent, Agent)
    assert isinstance(agent.activity, ActivityCondition)
    assert agent.activity.activity_type == 'activity'
    assert agent.activity.is_active


def test_phosphorylation_one_site_with_evidence():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    g = pybel.BELGraph()
    ev_text = 'Some evidence.'
    ev_pmid = '123456'
    edge_hash = g.add_qualified_edge(mek, erk, relation=pc.DIRECTLY_INCREASES,
                                     evidence=ev_text, citation=ev_pmid,
                                     annotations={"TextLocation": 'Abstract'})
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], Phosphorylation)
    assert pbp.statements[0].residue == 'T'
    assert pbp.statements[0].position == '185'
    enz = pbp.statements[0].enz
    sub = pbp.statements[0].sub
    assert enz.name == 'MAP2K1'
    assert enz.mods == []
    assert sub.name == 'MAPK1'
    assert sub.mods == []
    # Check evidence
    assert len(pbp.statements[0].evidence) == 1
    ev = pbp.statements[0].evidence[0]
    assert ev.source_api == 'pybel'
    assert ev.source_id == edge_hash
    assert ev.pmid == ev_pmid
    assert ev.text == ev_text
    assert ev.annotations == {'bel': 'p(HGNC:MAP2K1) directlyIncreases p(HGNC:MAPK1, pmod(Ph, Thr, 185))'}
    assert ev.epistemics == {'direct': True, 'section_type': 'abstract'}


def test_phosphorylation_two_sites():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr'),
                            pmod('Ph', position=187, code='Tyr')])
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.DIRECTLY_INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 2
    stmt1 = pbp.statements[0]
    stmt2 = pbp.statements[1]
    assert stmt1.residue == 'T'
    assert stmt1.position == '185'
    assert stmt2.residue == 'Y'
    assert stmt2.position == '187'
    assert stmt1.sub.mods == []
    assert stmt2.sub.mods == []
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount1_prot_obj():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount2_rna_obj():
    # FIXME: Create a transcription-specific statement for p->rna
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = rna(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount3_deg():
    # FIXME: Create a stability-specific statement for p->deg(p(Foo))
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         object_modifier=degradation(),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], DecreaseAmount)
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount4_subj_act():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         subject_modifier=activity(name='tscript'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert isinstance(subj.activity, ActivityCondition)
    assert subj.activity.activity_type == 'transcription'
    assert subj.activity.is_active == True
    assert len(pbp.statements[0].evidence) == 1

    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         subject_modifier=activity(name='act'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert isinstance(subj.activity, ActivityCondition)
    assert subj.activity.activity_type == 'activity'
    assert subj.activity.is_active == True
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_activity():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         subject_modifier=activity(name='kin'),
                         object_modifier=activity(name='kin'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], Activation)
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert isinstance(subj.activity, ActivityCondition)
    assert subj.activity.activity_type == 'kinase'
    assert subj.activity.is_active == True
    obj = pbp.statements[0].obj
    assert obj.name == 'MAPK1'
    assert obj.activity is None
    assert pbp.statements[0].obj_activity == 'kinase'
    assert len(pbp.statements[0].evidence) == 1


def test_active_form():
    p53_pmod = protein(name='TP53', namespace='HGNC',
                       variants=[pmod('Ph', position=33, code='Ser')])
    p53_obj = protein(name='TP53', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(p53_pmod, p53_obj, relation=pc.INCREASES,
                         object_modifier=activity(name='tscript'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, ActiveForm)
    assert stmt.activity == 'transcription'
    assert stmt.is_active is True
    ag = stmt.agent
    assert ag.name == 'TP53'
    assert len(ag.mods) == 1
    mc = ag.mods[0]
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'S'
    assert mc.position == '33'
    assert len(pbp.statements[0].evidence) == 1


def test_gef():
    sos = protein(name='SOS1', namespace='HGNC')
    kras = protein(name='KRAS', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(sos, kras, relation=pc.DIRECTLY_INCREASES,
                         subject_modifier=activity(name='activity'),
                         object_modifier=activity(name='gtp'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Gef)
    assert stmt.gef.name == 'SOS1'
    assert stmt.ras.name == 'KRAS'
    assert stmt.gef.activity.activity_type == 'activity'
    assert stmt.gef.activity.is_active is True
    assert stmt.ras.activity is None
    assert len(pbp.statements[0].evidence) == 1


def test_indirect_gef_is_activation():
    sos = protein(name='SOS1', namespace='HGNC')
    kras = protein(name='KRAS', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(sos, kras, relation=pc.INCREASES,
                         subject_modifier=activity(name='activity'),
                         object_modifier=activity(name='gtp'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Activation)
    assert stmt.subj.name == 'SOS1'
    assert stmt.obj.name == 'KRAS'
    assert stmt.subj.activity.activity_type == 'activity'
    assert stmt.subj.activity.is_active is True
    assert stmt.obj.activity is None
    assert stmt.obj_activity == 'gtpbound'
    assert len(pbp.statements[0].evidence) == 1


def test_gap():
    sos = protein(name='RASA1', namespace='HGNC')
    kras = protein(name='KRAS', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(sos, kras, relation=pc.DIRECTLY_DECREASES,
                         subject_modifier=activity(name='activity'),
                         object_modifier=activity(name='gtp'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Gap)
    assert stmt.gap.name == 'RASA1'
    assert stmt.ras.name == 'KRAS'
    assert stmt.gap.activity.activity_type == 'activity'
    assert stmt.gap.activity.is_active is True
    assert stmt.ras.activity is None
    assert len(pbp.statements[0].evidence) == 1


def test_activation_bioprocess():
    bax = protein(name='BAX', namespace='HGNC')
    apoptosis = bioprocess(name='apoptosis', namespace='GOBP')
    g = pybel.BELGraph()
    g.add_qualified_edge(bax, apoptosis, relation=pc.INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Activation)
    assert stmt.subj.name == 'BAX'
    assert stmt.obj.name == 'apoptosis'
    assert stmt.obj.db_refs == {} # FIXME: Update when GO lookup is implemented
    assert len(pbp.statements[0].evidence) == 1


def test_gtpactivation():
    kras = protein(name='KRAS', namespace='HGNC')
    braf = protein(name='BRAF', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(kras, braf, relation=pc.DIRECTLY_INCREASES,
                         subject_modifier=activity(name='gtp'),
                         object_modifier=activity(name='kin'),
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, GtpActivation)
    assert stmt.subj.name == 'KRAS'
    assert stmt.subj.activity.activity_type == 'gtpbound'
    assert stmt.subj.activity.is_active is True
    assert stmt.obj.name == 'BRAF'
    assert stmt.obj.activity is None
    assert stmt.obj_activity == 'kinase'
    assert len(stmt.evidence) == 1

if __name__ == '__main__':
    test_gtpactivation()
