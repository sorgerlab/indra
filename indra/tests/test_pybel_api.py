# -*- coding: utf-8 -*-

"""Tests for the PyBEL processor."""

import os
from urllib import request

from pybel import BELGraph
from pybel.dsl import *
from pybel.language import Entity
from pybel.io import from_nodelink_file
from pybel.examples import egf_graph
from indra.statements import *
from indra.sources import bel
from indra.sources.bel import processor as pb
from indra.sources.bel.api import process_cbn_jgif_file, process_pybel_graph, \
    small_corpus_url
from indra.databases import hgnc_client
from indra.statements.validate import assert_valid_statement

mek_hgnc_id = hgnc_client.get_hgnc_id('MAP2K1')
mek_up_id = hgnc_client.get_uniprot_id(mek_hgnc_id)


def test_pybel_neighborhood_query():
    bp = bel.process_pybel_neighborhood(['TP63'],
                                        network_type='graph_jsongz_url',
                                        network_file=small_corpus_url)
    assert bp.statements
    for stmt in bp.statements:
        assert_valid_statement(stmt)
    assert all([s.evidence[0].context is not None
                for s in bp.statements])
    assert all([s.evidence[0].context.cell_line.name == 'MCF 10A'
               for s in bp.statements])
    # Locate statement about epidermis development
    stmt = [st for st in bp.statements if st.agent_list()[1].name ==
            'epidermis development'][0]
    assert repr(stmt.evidence[0].context) == str(stmt.evidence[0].context)
    assert stmt.evidence[0].context == BioContext(
        location=RefContext(name="Cytoplasm",
                            db_refs={'MESH': 'D003593'}),
        cell_line=RefContext(name="MCF 10A",
                             db_refs={'EFO': '0001200'}),
        cell_type=RefContext(name="keratinocyte",
                             db_refs={'CL': 'CL:0000312'}),
        organ=RefContext(name="colon",
                         db_refs={'UBERON': 'UBERON:0001155'}),
        disease=RefContext(name="cancer",
                           db_refs={'DOID': 'DOID:162'}),
        species=RefContext(name="Rattus norvegicus",
                           db_refs={'TAXONOMY': '10116'})), \
        stmt.evidence[0].context
    # Test annotation manager
    assert bp.annot_manager.get_mapping('Species', '9606') == \
        'Homo sapiens'


def test_pybel_readme_example():
    bel_processor = bel.process_pybel_neighborhood(['KRAS', 'BRAF'])
    assert bel_processor.statements


def test_process_pybel():
    pbp = bel.process_pybel_graph(egf_graph)
    assert pbp.statements


def test_process_jgif():
    test_file_url = 'https://s3.amazonaws.com/bigmech/travis/Hox-2.0-Hs.jgf'
    test_file = 'Hox-2.0-Hs.jgf'
    if not os.path.exists(test_file):
        request.urlretrieve(url=test_file_url, filename=test_file)
    pbp = process_cbn_jgif_file(test_file)

    # Clean up
    os.remove(test_file)

    assert len(pbp.statements) == 26, len(pbp.statements)
    assert isinstance(pbp.statements[0], Statement)
    assert all(s.evidence[0].source_api == 'bel' for s in pbp.statements)


def test_nodelink_json():
    test_file_url = \
        'https://s3.amazonaws.com/bigmech/travis/Hox-2.0-Hs_nljson.json'
    test_file = 'Hox-2.0-Hs_nljson.json'
    if not os.path.exists(test_file):
        request.urlretrieve(url=test_file_url, filename=test_file)
    pbp = process_pybel_graph(from_nodelink_file(test_file))

    # Clean up
    os.remove(test_file)

    # Changed to 24, not really sure how to debug this one
    assert len(pbp.statements) == 24, (len(pbp.statements), pbp.statements)
    assert isinstance(pbp.statements[0], Statement)
    assert all(s.evidence[0].source_api == 'bel' for s in pbp.statements)


def test_get_agent_hgnc():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1', agent
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id

    # Now create an agent with an identifier
    mek = Protein(name='Foo', namespace='HGNC', identifier='6840')
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1', agent
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id


def test_get_agent_up():
    mek = Protein(namespace='UP', identifier='Q02750')
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == mek_hgnc_id
    assert agent.db_refs.get('UP') == mek_up_id


def test_get_agent_egid():
    node_data = Protein(name='5008', namespace='EGID')
    agent = pb.get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'OSM'
    assert len(agent.db_refs) == 3
    assert agent.db_refs['EGID'] == '5008'
    assert agent.db_refs['HGNC'] == '8506'
    assert agent.db_refs['UP'] == 'P13725'


def test_get_agent_mgi():
    node = Protein(namespace='MGI', name='Nr1h3')
    agent = pb.get_agent(node, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'Nr1h3'
    assert len(agent.db_refs) == 1
    assert agent.db_refs.get('UP') == 'Q9Z0Y9', agent.db_refs


def test_get_agent_rgd():
    node = Protein(namespace='RGD', name='Tp53')
    agent = pb.get_agent(node, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'Tp53'
    assert len(agent.db_refs) == 1
    assert agent.db_refs.get('UP') == 'P10361', agent.db_refs


def test_get_agent_sfam():
    node_data = Protein(
        namespace='SFAM',
        name='PRKC Family',
    )
    agent = pb.get_agent(node_data)
    assert isinstance(agent, Agent)
    assert len(agent.db_refs) == 2, agent.db_refs
    assert agent.db_refs['SFAM'] == 'F0198', agent.db_refs
    assert agent.db_refs['FPLX'] == 'PKC'
    assert agent.name == 'PKC'


def test_get_agent_sdis():
    node_data = Pathology(namespace='SDIS', name='metastasis')
    agent = pb.get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'metastasis'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['SDIS'] == 'D0340', agent.db_refs


def test_get_agent_chebi():
    node_data = Abundance(namespace='CHEBI', name='nitric oxide')
    agent = pb.get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'nitric oxide'
    assert len(agent.db_refs) == 1
    assert agent.db_refs['CHEBI'] == 'CHEBI:16480'


def test_get_agent_schem():
    node_data = Abundance(namespace='SCHEM', name='Promegestone')
    agent = pb.get_agent(node_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'Promegestone'
    assert len(agent.db_refs) == 2, agent.db_refs
    assert agent.db_refs['SCHEM'] == 'A2173', agent.db_refs


def test_get_agent_mirna():
    m = MicroRna(namespace='HGNC', name='MIRLET7A1')
    agent = pb.get_agent(m, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MIRLET7A1'
    assert agent.db_refs.get('MIRBASE') == 'MI0000060'
    assert agent.db_refs.get('HGNC') == '31476'

    m = MicroRna(namespace='HGNC', name='MIRLET7A1', identifier='31476')
    agent = pb.get_agent(m, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MIRLET7A1'
    assert agent.db_refs.get('MIRBASE') == 'MI0000060'
    assert agent.db_refs.get('HGNC') == '31476'

    m = MicroRna(namespace='MIRBASE', name='hsa-let-7a-1')
    agent = pb.get_agent(m, {})
    assert isinstance(agent, Agent)
    assert agent.name == 'MIRLET7A1'
    assert agent.db_refs.get('MIRBASE') == 'MI0000060'
    assert agent.db_refs.get('HGNC') == '31476'


def test_get_agent_fusion():
    node_data = ProteinFusion(
        partner_5p=Protein(namespace='HGNC', name='BCR'),
        partner_3p=Protein(namespace='HGNC', name='ABL1'),
    )
    agent = pb.get_agent(node_data)
    assert agent is None


def test_get_agent_up_no_id():
    mek = Protein(name='MAP2K1', namespace='UP')
    agent = pb.get_agent(mek, {})
    assert agent is None


def test_get_agent_meshpp():
    apoptosis = bioprocess(name='Apoptosis', namespace='MESHPP')
    agent = pb.get_agent(apoptosis)
    assert isinstance(agent, Agent)
    assert agent.name == 'Apoptosis'
    assert 'MESH' in agent.db_refs


def test_get_agent_meshd():
    hyperoxia = bioprocess(name='Hyperoxia', namespace='MESHD')
    agent = pb.get_agent(hyperoxia)
    assert isinstance(agent, Agent)
    assert agent.name == 'Hyperoxia'
    assert 'MESH' in agent.db_refs


def test_get_agent_with_mods():
    mek = Protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph')])
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert not mod.residue
    assert not mod.position

    mek = Protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', code='Ser')])
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert mod.residue == 'S'
    assert not mod.position

    mek = Protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', position=218)])
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert not mod.residue
    assert mod.position == '218'

    mek = Protein(name='MAP2K1', namespace='HGNC',
                  variants=[pmod('Ph', position=218, code='Ser')])
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mods) == 1
    mod = agent.mods[0]
    assert mod.mod_type == 'phosphorylation'
    assert mod.residue == 'S'
    assert mod.position == '218'


def test_get_agent_with_muts():
    mek = Protein(name='MAP2K1', namespace='HGNC',
                  variants=[hgvs('p.Val600Glu')])
    agent = pb.get_agent(mek, {})
    assert isinstance(agent, Agent)
    assert len(agent.mutations) == 1
    mut = agent.mutations[0]
    assert mut.position == '600'
    assert mut.residue_from == 'V'
    assert mut.residue_to == 'E'


def test_get_agent_with_activity():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    agent = pb.get_agent(mek, activity('act'))
    assert isinstance(agent, Agent)
    assert isinstance(agent.activity, ActivityCondition)
    assert agent.activity.activity_type == 'activity'
    assert agent.activity.is_active


def test_get_agent_complex():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb.get_agent(cplx)
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
    mek = Protein(name='MAP2K1', namespace='FOO')
    erk = Protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb.get_agent(cplx)
    assert agent is None

    # Bound agent is None
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='FOO',
                  variants=[pmod('Ph', position=185, code='Thr')])
    cplx = complex_abundance([mek, erk])
    agent = pb.get_agent(cplx)
    assert agent is None


def test_get_agent_named_complex_go():
    # TODO: Handle named complexes and map to FamPlex where possible
    node_data = NamedComplexAbundance(namespace='GOCCID', name='0043509')
    agent = pb.get_agent(node_data)
    assert agent is None


def test_get_agent_with_translocation():
    node_data = Protein(name='MAPK1', namespace='HGNC')
    # Some example edge data
    edge_data = translocation(
        from_loc=Entity(namespace='GOCC', name='intracellular'),
        to_loc=Entity(namespace='GOCC', name='extracellular space'),
    )
    agent = pb.get_agent(node_data, edge_data)
    assert isinstance(agent, Agent)
    assert agent.name == 'MAPK1'
    assert agent.location == 'extracellular space'


def test_phosphorylation_one_site_with_evidence():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr')])
    g = BELGraph()
    g.annotation_list['TextLocation'] = {'Abstract'}
    ev_text = 'Some evidence.'
    ev_pmid = '123456'
    edge_hash = g.add_directly_increases(
        mek, erk, evidence=ev_text,
        citation=ev_pmid,
        annotations={"TextLocation": 'Abstract'},
    )
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
    assert ev.pmid == ev_pmid, (ev.pmid, ev_pmid)
    assert ev.text == ev_text
    assert ev.annotations == {
        'bel': 'p(HGNC:MAP2K1) directlyIncreases '
               'p(HGNC:MAPK1, pmod(go:0006468 ! "protein phosphorylation", Thr, 185))'
    }
    assert ev.epistemics == {'direct': True, 'section_type': 'abstract'}


def test_doi_evidence():
    """Test processing edges with DOI citations."""
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.annotation_list['TextLocation'] = {'Abstract'}
    ev_doi = '123456'
    g.add_directly_increases(
        mek, erk, evidence='Some evidence.',
        citation=('doi', ev_doi),
        annotations={"TextLocation": 'Abstract'},
    )
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert len(pbp.statements[0].evidence) == 1
    ev = pbp.statements[0].evidence[0]
    assert ev.pmid is None
    assert 'DOI' in ev.text_refs
    assert ev.text_refs['DOI'] == ev_doi


def test_phosphorylation_two_sites():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC',
                  variants=[pmod('Ph', position=185, code='Thr'),
                            pmod('Ph', position=187, code='Tyr')])
    g = BELGraph()
    g.add_directly_increases(mek, erk, evidence="Some evidence.",
                             citation='123456')
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
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.add_increases(mek, erk, evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount2_rna_obj():
    # FIXME: Create a transcription-specific statement for p->rna
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = rna(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.add_increases(mek, erk, evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount3_deg():
    # FIXME: Create a stability-specific statement for p->deg(p(Foo))
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.add_increases(mek, erk, object_modifier=degradation(),
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], DecreaseAmount), pbp.statements[0]
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_amount4_subj_act():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.add_increases(mek, erk, source_modifier=activity(name='tscript'),
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert subj.activity is not None
    assert isinstance(subj.activity, ActivityCondition), subj.activity.__class__
    assert subj.activity.activity_type == 'transcription'
    assert subj.activity.is_active
    assert len(pbp.statements[0].evidence) == 1

    g = BELGraph()
    g.add_increases(mek, erk, source_modifier=activity(),
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert isinstance(subj.activity, ActivityCondition)
    assert subj.activity.activity_type == 'activity'
    assert subj.activity.is_active
    assert len(pbp.statements[0].evidence) == 1


def test_regulate_activity():
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    g.add_increases(mek, erk, source_modifier=activity(name='kin'),
                    target_modifier=activity(name='kin'),
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], Activation), pbp.statements[0].__class__
    subj = pbp.statements[0].subj
    assert subj.name == 'MAP2K1'
    assert isinstance(subj.activity, ActivityCondition)
    assert subj.activity.activity_type == 'kinase'
    assert subj.activity.is_active
    obj = pbp.statements[0].obj
    assert obj.name == 'MAPK1'
    assert obj.activity is None
    assert pbp.statements[0].obj_activity == 'kinase'
    assert len(pbp.statements[0].evidence) == 1


def test_active_form():
    p53_pmod = Protein(name='TP53', namespace='HGNC',
                       variants=[pmod('Ph', position=33, code='Ser')])
    p53_obj = Protein(name='TP53', namespace='HGNC')
    g = BELGraph()
    g.add_increases(p53_pmod, p53_obj,
                    target_modifier=activity(name='tscript'),
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
    sos = Protein(name='SOS1', namespace='HGNC')
    kras = Protein(name='KRAS', namespace='HGNC')
    g = BELGraph()
    g.add_directly_increases(sos, kras,
                             source_modifier=activity(),
                             target_modifier=activity(name='gtp'),
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
    sos = Protein(name='SOS1', namespace='HGNC')
    kras = Protein(name='KRAS', namespace='HGNC')
    g = BELGraph()
    g.add_increases(sos, kras, source_modifier=activity(),
                    target_modifier=activity(name='gtp'),
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
    sos = Protein(name='RASA1', namespace='HGNC')
    kras = Protein(name='KRAS', namespace='HGNC')
    g = BELGraph()
    g.add_directly_decreases(sos, kras,
                             source_modifier=activity(),
                             target_modifier=activity(name='gtp'),
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
    bax = Protein(name='BAX', namespace='HGNC')
    apoptosis = bioprocess(name='apoptotic process', namespace='GOBP')
    g = BELGraph()
    g.add_increases(bax, apoptosis, evidence="Some evidence.",
                    citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Activation)
    assert stmt.subj.name == 'BAX'
    assert stmt.obj.name == 'apoptotic process'
    assert 'GO' in stmt.obj.db_refs
    assert len(pbp.statements[0].evidence) == 1


def test_gtpactivation():
    kras = Protein(name='KRAS', namespace='HGNC')
    braf = Protein(name='BRAF', namespace='HGNC')
    g = BELGraph()
    g.add_directly_increases(kras, braf,
                             source_modifier=activity(name='gtp'),
                             target_modifier=activity(name='kin'),
                             evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, GtpActivation), stmt
    assert stmt.subj.name == 'KRAS'
    assert stmt.subj.activity.activity_type == 'gtpbound'
    assert stmt.subj.activity.is_active is True
    assert stmt.obj.name == 'BRAF'
    assert stmt.obj.activity is None
    assert stmt.obj_activity == 'kinase'
    assert len(stmt.evidence) == 1


def test_conversion():
    enz = Protein(name='PLCG1', namespace='HGNC')
    react_1 = abundance('SCHEM',
                        '1-Phosphatidyl-D-myo-inositol 4,5-bisphosphate')
    p1 = abundance('SCHEM', 'Diacylglycerol')
    p2 = abundance('SCHEM', 'Inositol 1,4,5-trisphosphate')

    rxn = reaction(
        reactants=react_1,
        products=[p1, p2],
    )
    g = BELGraph()
    g.add_directly_increases(enz, rxn,
                             source_modifier=activity(),
                             evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, Conversion)
    assert stmt.subj.name == 'PLCG1'
    assert stmt.subj.activity is not None
    assert stmt.subj.activity.activity_type is not None
    assert stmt.subj.activity.activity_type == 'activity', f'Got: {stmt.subj.activity.activity_type}'
    assert stmt.subj.activity.is_active is True
    assert len(stmt.obj_from) == 1
    assert isinstance(stmt.obj_from[0], Agent)
    assert stmt.obj_from[0].name == '1-Phosphatidyl-D-myo-inositol ' \
                                    '4,5-bisphosphate'
    assert len(stmt.obj_to) == 2
    # why do these not appear in alphabetical order?
    # PyBEL sorts the nodes based on their BEL, and
    # Inositol 1,4,5-trisphosphate gets quoted.
    assert stmt.obj_to[0].name == 'Inositol 1,4,5-trisphosphate'
    assert stmt.obj_to[1].name == 'Diacylglycerol'
    assert len(stmt.evidence) == 1


def test_controlled_transloc_loc_cond():
    """Controlled translocations are currently not handled."""
    subj = Protein(name='MAP2K1', namespace='HGNC')
    obj = Protein(name='MAPK1', namespace='HGNC')
    g = BELGraph()
    transloc = translocation(
        from_loc=Entity(namespace='GOCC', name='intracellular'),
        to_loc=Entity(namespace='GOCC', name='extracellular space'),
    )
    g.add_increases(subj, obj, target_modifier=transloc,
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert not pbp.statements, pbp.statements


def test_subject_transloc_loc_cond():
    """Translocations of the subject are treated as location conditions on the
    subject (using the to_loc location as the condition)"""
    subj = Protein(name='MAP2K1', namespace='HGNC')
    obj = Protein(name='MAPK1', namespace='HGNC')
    transloc = translocation(
        from_loc=Entity(namespace='GOCC', name='intracellular'),
        to_loc=Entity(namespace='GOCC', name='extracellular space'),
    )
    g = BELGraph()
    g.add_increases(subj, obj, source_modifier=transloc,
                    evidence="Some evidence.", citation='123456')
    pbp = bel.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    stmt = pbp.statements[0]
    assert isinstance(stmt, IncreaseAmount)
    assert stmt.subj.name == 'MAP2K1'
    assert stmt.subj.location is not None
    assert stmt.subj.location == 'extracellular space'
    assert stmt.obj.name == 'MAPK1'


def test_subject_transloc_active_form():
    """ActiveForms where the subject is a translocation--should draw on the
    to-location of the subject."""
    subj = Protein(name='MAP2K1', namespace='HGNC')
    obj = Protein(name='MAP2K1', namespace='HGNC')
    transloc = translocation(
        from_loc=Entity(namespace='GOCC', name='intracellular'),
        to_loc=Entity(namespace='GOCC', name='extracellular space'),
    )
    g = BELGraph()
    g.add_increases(subj, obj, source_modifier=transloc,
                    target_modifier=activity(name='kin'),
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
    raf = Protein(name='BRAF', namespace='HGNC')
    mek = Protein(name='MAP2K1', namespace='HGNC')
    erk = Protein(name='MAPK1', namespace='HGNC')
    cplx = complex_abundance([raf, mek])
    g = BELGraph()
    g.add_directly_increases(cplx, erk,
                             target_modifier=activity(name='kin'),
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


def test_process_bel_stmts():
    bp = bel.process_bel_stmt('p(HGNC:MDM2) directlyDecreases '
                              'tscript(p(HGNC:TP53))')
    assert len(bp.statements) == 1
    assert isinstance(bp.statements[0], Inhibition), bp.statements
    assert bp.statements[0].subj.name == 'MDM2', bp.statements
    assert bp.statements[0].obj.name == 'TP53', bp.statements

    bp = bel.process_bel_stmt('a(CHEBI:lipoprotein) increases '
                              'bp(GOBP:"inflammatory response")')
    assert len(bp.statements) == 1
    assert isinstance(bp.statements[0], Activation), bp.statements
    assert bp.statements[0].subj.name == 'lipoprotein', bp.statements
    assert bp.statements[0].obj.name == 'inflammatory response', bp.statements
