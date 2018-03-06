import pybel
from pybel.dsl import *
import pybel.constants as pc
from pybel.examples import egf_graph, sialic_acid_graph
from indra.statements import *
from indra.sources import bel
from indra.sources.bel import pybel_processor as pb
from indra.databases import hgnc_client


mek_hgnc_id = hgnc_client.get_hgnc_id('MAP2K1')
mek_up_id = hgnc_client.get_uniprot_id(mek_hgnc_id)


def test_process_pybel():
    pbp = bel.process_pybel_graph(egf_graph)
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


def test_get_agent_egid():
    node_data = {'function': 'Protein', 'name': '5008', 'namespace': 'EGID'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'OSM'
    assert len(agent.db_refs) == 3
    assert agent.db_refs['EGID'] == '5008'
    assert agent.db_refs['HGNC'] == '8506'
    assert agent.db_refs['UP'] == 'P13725'


def test_get_agent_mgi():
    node = protein(namespace='MGI', name='Nr1h3')
    agent = pb._get_agent(node, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'Nr1h3'
    assert len(agent.db_refs) == 1
    assert agent.db_refs.get('MGI') == 'Nr1h3'


def test_get_agent_rgd():
    node = protein(namespace='RGD', name='Tp53')
    agent = pb._get_agent(node, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'Tp53'
    assert len(agent.db_refs) == 1
    assert agent.db_refs.get('RGD') == 'Tp53'


def test_get_agent_sfam():
    node_data = {
            'cname': 'PRKC Family',
            'function': 'Protein',
            'name': 'PRKC Family',
            'namespace': 'SFAM'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert len(agent.db_refs) == 2
    assert agent.db_refs['SFAM'] == 'PRKC Family'
    assert agent.db_refs['FPLX'] == 'PKC'
    assert agent.name == 'PKC'


def test_get_agent_sdis():
    node_data = {
            'cname': 'metastasis',
            'function': 'Pathology',
            'name': 'metastasis',
            'namespace': 'SDIS'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'metastasis'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['SDIS'] == 'metastasis'


def test_get_agent_chebi():
    node_data = {
            'cname': 'nitric oxide',
            'function': 'Abundance',
            'name': 'nitric oxide',
            'namespace': 'CHEBI'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'nitric oxide'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['CHEBI'] == 'CHEBI:16480'


def test_get_agent_schem():
    node_data = {
            'cname': 'Promegestone',
            'function': 'Abundance',
            'name': 'Promegestone',
            'namespace': 'SCHEM'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'Promegestone'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['SCHEM'] == 'Promegestone'


def test_get_agent_mirna():
    node_data = {
            'cname': 'MIR218-1',
            'function': 'miRNA',
            'name': 'MIR218-1',
            'namespace': 'HGNC'}
    agent = pb._get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'MIR218-1'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['HGNC'] == '31595'


def test_get_agent_fusion():
    node_data = {'function': 'Protein',
                 'fusion': {
                     'partner_5p': {'namespace': 'HGNC', 'name': 'BCR'},
                     'range_5p': {'missing': '?'},
                     'range_3p': {'missing': '?'},
                     'partner_3p': {'namespace': 'HGNC', 'name': 'ABL1'}}}
    agent = pb._get_agent(node_data)
    assert agent is None


def test_get_agent_up_no_id():
    mek = protein(name='MAP2K1', namespace='UP')
    agent = pb._get_agent(mek, {})
    assert agent is None


def test_get_agent_meshpp():
    apoptosis = bioprocess(name='Apoptosis', namespace='MESHPP')
    agent = pb._get_agent(apoptosis)
    assert isinstance(agent, Agent)
    assert agent.name == 'Apoptosis'
    assert agent.db_refs == {}


def test_get_agent_meshd():
    hyperoxia = bioprocess(name='Hyperoxia', namespace='MESHD')
    agent = pb._get_agent(hyperoxia)
    assert isinstance(agent, Agent)
    assert agent.name == 'Hyperoxia'
    assert agent.db_refs == {}


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


def test_get_agent_complex():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb._get_agent(cplx)
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert len(agent.bound_conditions) == 1
    bc = agent.bound_conditions[0]
    assert isinstance(bc, BoundCondition)
    assert bc.is_bound is True
    bc_agent = bc.agent
    assert bc_agent.name == 'MAPK1'
    assert len(bc_agent.mods) == 1
    assert bc_agent.mods[0].mod_type == 'phosphorylation'
    assert bc_agent.mods[0].residue == 'T'
    assert bc_agent.mods[0].position == '185'


def test_get_agent_complex_none_agent():
    """If one of the agents in the complex can't be obtained (e.g., an
    unhandled namespace), then the complex itself should be None."""
    # Prime agent is None
    mek = protein(name='MAP2K1', namespace='FOO')
    erk = protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb._get_agent(cplx)
    assert agent is None

    # Bound agent is None
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='FOO',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb._get_agent(cplx)
    assert agent is None


def test_get_agent_named_complex_go():
    # TODO: Handle named complexes and map to FamPlex where possible
    node_data = {
            'cname': '0043509',
            'function': 'Complex',
            'name': '0043509',
            'namespace': 'GOCCID'}
    agent = pb._get_agent(node_data)
    assert agent is None


def test_get_agent_with_translocation():
    node_data = protein(name='MAPK1', namespace='HGNC')
    # Some example edge data
    edge_data = translocation(from_loc=entity('GOCC', 'intracellular'),
                             to_loc=entity('GOCC', 'extracellular space'))
    agent = pb._get_agent(node_data, edge_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'MAPK1'
    assert agent.location == 'extracellular space'


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
    pbp = bel.process_pybel_graph(g)
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
    assert ev.source_api == 'bel'
    assert ev.source_id == edge_hash
    assert ev.pmid == ev_pmid
    assert ev.text == ev_text
    assert ev.annotations == {'bel': 'p(HGNC:MAP2K1) directlyIncreases '
                                     'p(HGNC:MAPK1, pmod(Ph, Thr, 185))'}
    assert ev.epistemics == {'direct': True, 'section_type': 'abstract'}


def test_phosphorylation_two_sites():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr'),
                            pmod('Ph', position=187, code='Tyr')])
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.DIRECTLY_INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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
    pbp = bel.process_pybel_graph(g)
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


def test_conversion():
    enz = protein(name='PLCG1', namespace='HGNC')
    rxn = {
            'function': 'Reaction',
            'reactants': [
                {'function': 'Abundance',
                 'name': '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate',
                 'namespace': 'SCHEM'}
            ],
            'products': [
                {'function': 'Abundance',
                 'name': 'Diacylglycerol',
                 'namespace': 'SCHEM'},
                {'cname': 'Inositol 1,4,5-trisphosphate',
                 'function': 'Abundance',
                 'name': 'Inositol 1,4,5-trisphosphate',
                 'namespace': 'SCHEM'}
            ]
          }
    g = pybel.BELGraph()
    g.add_qualified_edge(enz, rxn, relation=pc.DIRECTLY_INCREASES,
                         subject_modifier=activity(name='activity'),
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Conversion)
    assert stmt.subj.name == 'PLCG1'
    assert stmt.subj.activity.activity_type == 'activity'
    assert stmt.subj.activity.is_active is True
    assert len(stmt.obj_from) == 1
    assert isinstance(stmt.obj_from[0], Agent)
    assert stmt.obj_from[0].name == '1-Phosphatidyl-D-myo-inositol ' \
                                    '4,5-bisphosphate'
    assert len(stmt.obj_to) == 2
    assert stmt.obj_to[0].name == 'Diacylglycerol'
    assert stmt.obj_to[1].name == 'Inositol 1,4,5-trisphosphate'
    assert len(stmt.evidence) == 1


def test_controlled_transloc_loc_cond():
    """Controlled translocations are currently not handled."""
    subj = protein(name='MAP2K1', namespace='HGNC')
    obj = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    transloc = translocation(from_loc=entity('GOCC', 'intracellular'),
                             to_loc=entity('GOCC', 'extracellular space'))
    g.add_qualified_edge(subj, obj, relation=pc.INCREASES,
                         object_modifier=transloc,
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert not pbp.statements


def test_subject_transloc_loc_cond():
    """Translocations of the subject are treated as location conditions on the
    subject (using the to_loc location as the condition)"""
    subj = protein(name='MAP2K1', namespace='HGNC')
    obj = protein(name='MAPK1', namespace='HGNC')
    transloc = translocation(from_loc=entity('GOCC', 'intracellular'),
                             to_loc=entity('GOCC', 'extracellular space'))
    g = pybel.BELGraph()
    g.add_qualified_edge(subj, obj, relation=pc.INCREASES,
                         subject_modifier=transloc,
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, IncreaseAmount)
    assert stmt.subj.name == 'MAP2K1'
    assert stmt.subj.location == 'extracellular space'
    assert stmt.obj.name == 'MAPK1'


def test_subject_transloc_active_form():
    """ActiveForms where the subject is a translocation--should draw on the
    to-location of the subject."""
    subj = protein(name='MAP2K1', namespace='HGNC')
    obj = protein(name='MAP2K1', namespace='HGNC')
    transloc = translocation(from_loc=entity('GOCC', 'intracellular'),
                             to_loc=entity('GOCC', 'extracellular space'))
    g = pybel.BELGraph()
    g.add_qualified_edge(subj, obj, relation=pc.INCREASES,
                         subject_modifier=transloc,
                         object_modifier=activity(name='kin'),
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, ActiveForm)
    assert stmt.agent.name == 'MAP2K1'
    assert stmt.agent.location == 'extracellular space'
    assert stmt.agent.activity is None
    assert stmt.activity == 'kinase'
    assert stmt.is_active is True


def test_complex_stmt_with_activation():
    raf = protein(name='BRAF', namespace='HGNC')
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    cplx = complex_abundance([raf, mek])
    g = pybel.BELGraph()
    g.add_qualified_edge(cplx, erk, relation=pc.DIRECTLY_INCREASES,
                         object_modifier=activity(name='kin'),
                         evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 2
    stmt1 = pbp.statements[0]
    assert isinstance(stmt1, Complex)
    assert len(stmt1.agent_list()) == 2
    assert sorted([ag.name for ag in stmt1.agent_list()]) == ['BRAF', 'MAP2K1']
    assert stmt1.evidence
    stmt2 = pbp.statements[1]
    assert isinstance(stmt2, Activation)
    assert stmt2.subj.name == 'BRAF'
    assert stmt2.subj.bound_conditions[0].agent.name == 'MAP2K1'
    assert stmt2.obj.name == 'MAPK1'
    assert stmt2.obj.activity is None
    assert stmt2.obj_activity == 'kinase'


if __name__ == '__main__':
    test_get_agent_fusion()

