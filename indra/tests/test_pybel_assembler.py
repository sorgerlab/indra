# -*- coding: utf-8 -*-

"""Tests for the PyBEL assembler."""

import json

import networkx as nx
import pybel.constants as pc
from pybel.dsl import abundance, activity, bioprocess, \
    complex_abundance, hgvs, pmod, protein, reaction

from indra.assemblers.pybel import assembler as pa
from indra.databases import hgnc_client
from indra.statements import *


def id(gene_name):
    return hgnc_client.get_hgnc_id(gene_name)


phos_dsl = pmod('Ph', 'Ser', 218)
ub_dsl = pmod('Ub', 'Ser', 218)
egfr_phos_dsl = pmod('Ph', 'Tyr', 1173)

braf_dsl = protein(namespace='HGNC', name='BRAF', identifier='1097')
map2k1_dsl = protein(namespace='HGNC', name='MAP2K1', identifier='6840')
tp53_dsl = protein(namespace='HGNC', name='TP53', identifier='11998')
mdm2_dsl = protein(namespace='HGNC', name='MDM2', identifier='6973')
egfr_dsl = protein(namespace='HGNC', name='EGFR', identifier='3236')

chebi_17534 = abundance(namespace='CHEBI', name='D-glucose',
                        identifier='17634')
chebi_4170 = abundance(namespace='CHEBI', name='D-glucopyranose 6-phosphate',
                       identifier='4170')
chebi_17534_to_4170 = reaction(chebi_17534, chebi_4170)

grb2_dsl = protein(namespace='HGNC', name='GRB2', identifier='4566')
sos1_dsl = protein(namespace='HGNC', name='SOS1', identifier='11187')
sos1_phosphorylated_dsl = sos1_dsl.with_variants(pmod('Ph'))
kras_node = protein(namespace='HGNC', name='KRAS', identifier='6407')

egfr_grb2_sos1_complex_dsl = complex_abundance([
    egfr_dsl,
    grb2_dsl,
    sos1_dsl,
])

egfr_grb2_sos1_phos_complex_dsl = complex_abundance([
    egfr_dsl,
    grb2_dsl,
    sos1_phosphorylated_dsl,
])


def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')


def get_edge_data(g, u, v):
    assert g.has_edge(u, v)
    data = g.get_edge_data(u, v)
    return list(data.values())[0]


def get_first_edge_data(g):
    return list(g.edges(data=True))[0][2]


def test_simple_modification_no_evidence():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_cat = Agent('BRAF', activity=ActivityCondition('catalytic', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})  # MEK
    stmt1 = Phosphorylation(braf, map2k1, 'S', '218')
    stmt2 = Phosphorylation(braf_kin, map2k1, 'S', '218')
    stmt3 = Ubiquitination(braf_cat, map2k1, 'S', '218')
    # Edge info for subject
    edge1 = None
    edge2 = activity('kin')
    edge3 = activity('cat')
    for stmt, modtuple, subj_edge in ((stmt1, phos_dsl, edge1),
                                      (stmt2, phos_dsl, edge2),
                                      (stmt3, ub_dsl, edge3)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert belgraph.number_of_nodes() == 3, belgraph.number_of_nodes()
        map2k1_mod_dsl = map2k1_dsl.with_variants(modtuple)
        assert set(belgraph) == {braf_dsl, map2k1_dsl, map2k1_mod_dsl}, \
            (set(belgraph), {braf_dsl, map2k1_dsl, map2k1_mod_dsl})
        assert belgraph.number_of_edges() == 2, belgraph.number_of_edges()
        assert belgraph.has_edge(map2k1_dsl, map2k1_mod_dsl)
        assert belgraph.has_edge(braf_dsl, map2k1_mod_dsl)
        edge_data = get_edge_data(belgraph, braf_dsl, map2k1_mod_dsl)
        assert edge_data[pc.RELATION] == pc.INCREASES
        assert edge_data.get(pc.SUBJECT) == subj_edge


def test_modification_with_evidences():
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    evidence = Evidence(source_api='test', text='evidence text', pmid='1234', epistemics={
        'dummy': ['a', 'b'],
        'scalar': 'yes',
        'missing': None,
    })
    stmt = Phosphorylation(braf_kin, mek, 'S', '218', evidence=evidence)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 3, belgraph.number_of_nodes()
    assert braf_dsl in belgraph
    map2k1_mod_dsl = map2k1_dsl.with_variants(phos_dsl)
    assert map2k1_mod_dsl in belgraph
    assert belgraph.number_of_edges() == 2
    edge_data = get_edge_data(belgraph, braf_dsl, map2k1_mod_dsl)
    assert edge_data.get(pc.SUBJECT) == activity('kin')
    assert edge_data[pc.RELATION] == pc.INCREASES
    assert edge_data.get(pc.EVIDENCE) == 'evidence text', edge_data
    assert edge_data[pc.CITATION] == {
        pc.CITATION_DB: pc.CITATION_TYPE_PUBMED,
        pc.CITATION_IDENTIFIER: '1234',
    }
    assert 'source_api' in edge_data[pc.ANNOTATIONS]
    assert 'test' in edge_data[pc.ANNOTATIONS]['source_api']
    assert 'source_id' not in edge_data[pc.ANNOTATIONS]
    assert 'source_hash' in edge_data[pc.ANNOTATIONS]
    assert 'dummy' in edge_data[pc.ANNOTATIONS]
    assert 'a' in edge_data[pc.ANNOTATIONS]['dummy']
    assert 'b' in edge_data[pc.ANNOTATIONS]['dummy']
    assert 'scalar' in edge_data[pc.ANNOTATIONS]
    assert 'yes' in edge_data[pc.ANNOTATIONS]['scalar']
    assert 'missing' not in edge_data[pc.ANNOTATIONS]


def test_modification_with_mutation():
    braf = Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')],
                 db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # Adds in the base protein nodes as well as the variants (so 4 nodes)
    assert belgraph.number_of_nodes() == 4, belgraph.number_of_nodes()
    braf_mut_dsl = braf_dsl.with_variants(hgvs('p.Val600Glu'))
    assert braf_mut_dsl in belgraph


def test_activation():
    braf_no_act = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt1 = Activation(braf_no_act, mek)
    stmt2 = Activation(braf_kin, mek, 'kinase')
    hash1 = stmt1.get_hash(refresh=True)
    hash2 = stmt2.get_hash(refresh=True)
    edge1 = {
        pc.RELATION: pc.INCREASES,
        pc.OBJECT: activity(),
        pc.ANNOTATIONS: {
            'stmt_hash': {hash1: True},
            'uuid': {stmt1.uuid: True},
            'belief': {stmt1.belief: True},
        },
    }
    edge2 = {
        pc.RELATION: pc.INCREASES,
        pc.SUBJECT: activity('kin'),
        pc.OBJECT: activity('kin'),
        pc.ANNOTATIONS: {
            'stmt_hash': {hash2: True},
            'uuid': {stmt2.uuid: True},
            'belief': {stmt2.belief: True},
        },
    }
    for stmt, edge in ((stmt1, edge1), (stmt2, edge2)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert belgraph.number_of_nodes() == 2, belgraph.number_of_nodes()
        assert braf_dsl in belgraph
        assert map2k1_dsl in belgraph
        assert belgraph.number_of_edges() == 1
        edge_data = get_first_edge_data(belgraph)
        assert edge_data == edge, edge_data


def test_direct_activation():
    braf_no_act = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt1_ev = Evidence(
        pmid='1234',
        epistemics={'direct': True},
    )
    stmt1 = Activation(braf_no_act, mek, evidence=stmt1_ev)
    stmt2 = Activation(braf_kin, mek, 'kinase', evidence=stmt1_ev)
    hash1 = stmt1.get_hash(refresh=True)
    hash2 = stmt2.get_hash(refresh=True)
    edge1 = {
        pc.RELATION: pc.DIRECTLY_INCREASES,
        pc.OBJECT: activity(),
        pc.EVIDENCE: 'No evidence text.',
        pc.CITATION: {
            pc.CITATION_DB: pc.CITATION_TYPE_PUBMED,
            pc.CITATION_IDENTIFIER: '1234',
        },
        pc.ANNOTATIONS: {
            'stmt_hash': {hash1: True},
            'source_hash': {stmt1_ev.get_source_hash(): True},
            'uuid': {stmt1.uuid: True},
            'belief': {stmt1.belief: True},
        },
    }
    edge2 = {
        pc.RELATION: pc.DIRECTLY_INCREASES,
        pc.SUBJECT: activity('kin'),
        pc.OBJECT: activity('kin'),
        pc.EVIDENCE: 'No evidence text.',
        pc.CITATION: {
            pc.CITATION_DB: pc.CITATION_TYPE_PUBMED,
            pc.CITATION_IDENTIFIER: '1234',
        },
        pc.ANNOTATIONS: {
            'stmt_hash': {hash2: True},
            'source_hash': {stmt1_ev.get_source_hash(): True},
            'uuid': {stmt2.uuid: True},
            'belief': {stmt2.belief: True},
        },
    }
    for stmt, expected_edge in ((stmt1, edge1), (stmt2, edge2)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert belgraph.number_of_nodes() == 2, belgraph.number_of_nodes()
        assert braf_dsl in belgraph
        assert map2k1_dsl in belgraph
        assert belgraph.number_of_edges() == 1
        edge_data = get_first_edge_data(belgraph)
        assert expected_edge == edge_data, json.dumps(edge_data, indent=1)


def test_inhibition():
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Inhibition(braf_kin, mek, 'kinase')
    stmt_hash = stmt.get_hash(refresh=True)
    edge = {
        pc.RELATION: pc.DECREASES,
        pc.SUBJECT: activity('kin'),
        pc.OBJECT: activity('kin'),
        pc.ANNOTATIONS: {
            'stmt_hash': {stmt_hash: True},
            'uuid': {stmt.uuid: True},
            'belief': {stmt.belief: True},
        },
    }
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 2, belgraph.number_of_nodes()
    assert braf_dsl in belgraph
    assert map2k1_dsl in belgraph
    assert belgraph.number_of_edges() == 1
    edge_data = get_first_edge_data(belgraph)
    assert edge_data == edge, edge_data


def test_increase_amount():
    tp53 = Agent('TP53', db_refs={'HGNC': '11998'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '6973'})

    stmt = IncreaseAmount(tp53, mdm2)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 2, belgraph.number_of_nodes()
    assert mdm2_dsl in belgraph
    assert tp53_dsl in belgraph
    assert belgraph.number_of_edges() == 1
    edge_data = get_first_edge_data(belgraph)
    assert edge_data[pc.RELATION] == pc.INCREASES


def test_increase_amount_tscript():
    tp53 = Agent('TP53', activity=ActivityCondition('transcription', True),
                 db_refs={'HGNC': '11998'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '6973'})

    stmt = IncreaseAmount(tp53, mdm2)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 2, belgraph.number_of_nodes()
    assert mdm2_dsl in belgraph
    assert tp53_dsl in belgraph
    assert belgraph.number_of_edges() == 1
    edge_data = get_first_edge_data(belgraph)
    assert edge_data[pc.RELATION] == pc.INCREASES
    assert edge_data[pc.SUBJECT] == activity('tscript')


def test_gef():
    gef = Agent('SOS1', mods=[ModCondition('phosphorylation')],
                db_refs={'HGNC': '11187'})
    ras = Agent('KRAS', db_refs={'HGNC': '6407'})
    stmt = Gef(gef, ras)
    stmt_hash = stmt.get_hash(refresh=True)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 3
    assert belgraph.number_of_edges() == 2

    gef_reference_node = protein(
        namespace='HGNC', name='SOS1', identifier='11187')
    gef_node = gef_reference_node.with_variants(pmod('Ph'))
    assert gef_reference_node in belgraph
    assert gef_node in belgraph
    assert kras_node in belgraph

    edge_data = get_edge_data(belgraph, gef_node, kras_node)
    edge = {
        pc.RELATION: pc.DIRECTLY_INCREASES,
        pc.SUBJECT: activity('gef'),
        pc.OBJECT: activity('gtp'),
        pc.ANNOTATIONS: {
            'stmt_hash': {stmt_hash: True},
            'uuid': {stmt.uuid: True},
            'belief': {stmt.belief: True},
        },
    }
    assert edge_data == edge, edge_data


def test_gap():
    gap = Agent('RASA1', mods=[ModCondition('phosphorylation')],
                db_refs={'HGNC': '9871'})
    ras = Agent('KRAS', db_refs={'HGNC': '6407'})
    stmt = Gap(gap, ras)
    stmt_hash = stmt.get_hash(refresh=True)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 3
    assert belgraph.number_of_edges() == 2

    gap_reference_node = protein(
        namespace='HGNC', name='RASA1', identifier='9871')
    gap_node = gap_reference_node.with_variants(pmod('Ph'))
    ras_node = protein(namespace='HGNC', name='KRAS', identifier='6407')

    assert gap_reference_node in belgraph
    assert gap_node in belgraph
    assert ras_node in belgraph
    edge_data = get_edge_data(belgraph, gap_node, ras_node)
    edge = {
        pc.RELATION: pc.DIRECTLY_DECREASES,
        pc.SUBJECT: activity('gap'),
        pc.OBJECT: activity('gtp'),
        pc.ANNOTATIONS: {
            'stmt_hash': {stmt_hash: True},
            'uuid': {stmt.uuid: True},
            'belief': {stmt.belief: True},
        },
    }
    assert edge_data == edge, edge_data


def test_active_form():
    ras = Agent('KRAS', mutations=[MutCondition('12', 'G', 'V')],
                db_refs={'HGNC': '6407'})
    mapk1_p = Agent('MAP2K1',
                    mods=[ModCondition('phosphorylation', 'T', '185')],
                    db_refs={'HGNC': hgnc_client.get_hgnc_id('MAP2K1')})
    mapk1_pp = Agent('MAP2K1',
                     mods=[ModCondition('phosphorylation', 'T', '185'),
                           ModCondition('phosphorylation', 'Y', '187')],
                     db_refs={'HGNC': hgnc_client.get_hgnc_id('MAP2K1')})
    stmt1 = ActiveForm(ras, 'gtpbound', True)
    stmt2 = ActiveForm(mapk1_p, 'kinase', True)
    stmt3 = ActiveForm(mapk1_pp, 'kinase', True)
    for i, stmt in enumerate((stmt1, stmt2, stmt3)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        if i == 2:
            assert len(belgraph) == 3, len(belgraph)
        else:
            assert len(belgraph) == 2, len(belgraph)


def test_complex():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    grb2 = Agent('GRB2', db_refs={'HGNC': id('GRB2')})
    sos = Agent('SOS1', db_refs={'HGNC': id('SOS1')})
    stmt = Complex([egfr, grb2, sos])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # The graph should contain the node for the complex as well as nodes
    # for all of the members
    assert len(belgraph) == 4
    assert egfr_grb2_sos1_complex_dsl in belgraph
    for member in egfr_grb2_sos1_complex_dsl.members:
        assert member in belgraph


def test_rxn_no_controller():
    glu = Agent('D-glucose', db_refs={'CHEBI': 'CHEBI:17634'})
    g6p = Agent('D-glucopyranose 6-phosphate', db_refs={'CHEBI': 'CHEBI:4170'})
    stmt = Conversion(None, [glu], [g6p])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # The graph should contain the node for the reaction as well as nodes
    # for all of the members
    assert chebi_17534_to_4170 in belgraph
    for reactant in chebi_17534_to_4170.reactants:
        assert reactant in belgraph

    assert belgraph.number_of_nodes() == 3, belgraph.number_of_nodes()
    # TODO check edge chebi_17534_to_4170 hasReactant chebi_17534

    for product in chebi_17534_to_4170.products:
        assert product in belgraph
    # TODO check edge chebi_17534_to_4170 hasProduct chebi_4170


def test_rxn_with_controller():
    hk1 = Agent('HK1', db_refs={'HGNC': id('HK1')})
    glu = Agent('D-glucose', db_refs={'CHEBI': 'CHEBI:17634'})
    g6p = Agent('D-glucopyranose 6-phosphate', db_refs={'CHEBI': 'CHEBI:4170'})
    stmt = Conversion(hk1, [glu], [g6p])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()

    # check the catalyst makes it
    assert protein(namespace='HGNC', name='HK1', identifier=id('HK1')) \
        in belgraph

    # The reaction data should be the same as before
    assert chebi_17534 in belgraph
    assert chebi_4170 in belgraph
    assert chebi_17534_to_4170 in belgraph

    # The graph should contain the node for the reaction as well as nodes
    # for all of the members
    assert belgraph.number_of_nodes() == 4, belgraph.number_of_nodes()


def test_autophosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    stmt = Autophosphorylation(egfr, 'Y', '1173')
    stmt_hash = stmt.get_hash(refresh=True)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 2
    assert egfr_dsl in belgraph
    egfr_phos_node = egfr_dsl.with_variants(egfr_phos_dsl)
    assert egfr_dsl in belgraph
    assert egfr_phos_node in belgraph
    assert belgraph.number_of_nodes() == 2
    assert belgraph.number_of_edges() == 2
    # There will be two edges between these nodes
    edge_dicts = list(belgraph.get_edge_data(egfr_dsl,
                                             egfr_phos_node).values())
    assert {pc.RELATION: pc.DIRECTLY_INCREASES,
            pc.ANNOTATIONS: {
                'stmt_hash': {stmt_hash: True},
                'uuid': {stmt.uuid: True},
                'belief': {stmt.belief: True},
            }} in edge_dicts

    # Test an autophosphorylation with a bound condition
    tab1 = Agent('TAB1', db_refs={'HGNC': id('TAB1')})
    p38_tab1 = Agent('MAPK14', bound_conditions=[BoundCondition(tab1)],
                     db_refs={'HGNC': id('MAPK14')})
    stmt = Autophosphorylation(p38_tab1, 'Y', '100')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 4
    assert belgraph.number_of_edges() == 4


def test_bound_condition():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    grb2 = Agent('GRB2', db_refs={'HGNC': id('GRB2')})
    ras = Agent('KRAS', db_refs={'HGNC': '6407'})
    sos1_bound = Agent(
        'SOS1', mods=[ModCondition('phosphorylation')],
        bound_conditions=[BoundCondition(egfr), BoundCondition(grb2)],
        db_refs={'HGNC': id('SOS1')}
    )
    stmt = Gef(sos1_bound, ras)
    stmt_hash = stmt.get_hash(refresh=True)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 6
    assert belgraph.number_of_edges() == 5
    # Don't bother to check the tuple, which is now generated by
    # PyBEL directly, but check the node data

    assert egfr_grb2_sos1_phos_complex_dsl in belgraph
    assert kras_node in belgraph
    assert (egfr_grb2_sos1_phos_complex_dsl, kras_node) in belgraph.edges()

    edge_data = (
        egfr_grb2_sos1_phos_complex_dsl,
        kras_node,
        {
             pc.RELATION: pc.DIRECTLY_INCREASES,
             pc.OBJECT: activity('gtp'),
             pc.ANNOTATIONS: {
                 'stmt_hash': {stmt_hash: True},
                 'uuid': {stmt.uuid: True},
                 'belief': {stmt.belief: True},
             },
        },
    )
    belgraph_edges = belgraph.edges(data=True)
    assert edge_data in belgraph_edges, belgraph_edges


def test_transphosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    egfr_dimer = Agent('EGFR', bound_conditions=[BoundCondition(egfr)],
                       db_refs={'HGNC': id('EGFR')})
    stmt = Transphosphorylation(egfr_dimer, 'Y', '1173')
    stmt_hash = stmt.get_hash(refresh=True)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 3
    assert belgraph.number_of_edges() == 3

    egfr_dimer_node = complex_abundance([egfr_dsl, egfr_dsl])
    egfr_phos_node = egfr_dsl.with_variants(pmod('Ph', 'Tyr', 1173))
    edge_data = get_edge_data(belgraph, egfr_dimer_node, egfr_phos_node)
    assert edge_data == {
        pc.RELATION: pc.DIRECTLY_INCREASES,
        pc.ANNOTATIONS: {
            'stmt_hash': {stmt_hash: True},
            'uuid': {stmt.uuid: True},
            'belief': {stmt.belief: True},
        },
    }, edge_data


"""
def test_translocation():
    foxo = Agent('FOXO1', db_refs={'HGNC': id('FOXO1')})
    stmt = Translocation(foxo, 'cytoplasm', 'nucleus')
    nuc_go = 'GO:0005634'
    cyto_go = 'GO:0005737'
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 1
"""


def test_complex_with_pmod():
    sos1_phos = Agent('SOS1',
                      mods=[ModCondition('phosphorylation', 'Y', '100')],
                      db_refs={'HGNC': id('SOS1')})
    grb2 = Agent('GRB2', db_refs={'HGNC': id('GRB2')})
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    stmt = Complex([sos1_phos, grb2, egfr])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert belgraph.number_of_nodes() == 5
    assert belgraph.number_of_edges() == 4

    egfr_grb2_sos_phos_tyr_100 = complex_abundance([
        egfr_dsl,
        grb2_dsl,
        sos1_dsl.with_variants(pmod('Ph', 'Tyr', 100))
    ])

    assert sos1_dsl in belgraph
    assert egfr_grb2_sos_phos_tyr_100 in belgraph
    for member in egfr_grb2_sos_phos_tyr_100.members:
        assert member in belgraph


def test_complex_with_complex():
    grb2 = Agent('GRB2', db_refs={'HGNC': id('GRB2')})
    egfr_grb2 = Agent('EGFR', db_refs={'HGNC': id('EGFR')},
                      bound_conditions=[BoundCondition(grb2)])
    sos1_phos = Agent('SOS1',
                      mods=[ModCondition('phosphorylation', 'Y', '100')],
                      db_refs={'HGNC': id('SOS1')})
    stmt = Complex([sos1_phos, egfr_grb2])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 6
    assert belgraph.number_of_edges() == 5

    egfr_grb2_complex = complex_abundance([egfr_dsl, grb2_dsl])
    egfr_grb2_complex_sos1_phos_complex = complex_abundance([
        egfr_grb2_complex,
        sos1_dsl.with_variants(pmod('Ph', 'Tyr', 100))
    ])

    assert egfr_grb2_complex in belgraph
    for member in egfr_grb2_complex.members:
        assert member in belgraph

    assert egfr_grb2_complex_sos1_phos_complex in belgraph
    for member in egfr_grb2_complex_sos1_phos_complex.members:
        assert member in belgraph


def test_no_activity_on_bioprocess():
    yfg_agent = Agent('PPP1R13L', db_refs={'HGNC': id('PPP1R13L')})
    apoptosis_agent = Agent('apoptotic process', db_refs={'GO': 'GO:0006915'})

    stmt = Activation(yfg_agent, apoptosis_agent)
    pba = pa.PybelAssembler([stmt])

    belgraph = pba.make_model()
    assert len(belgraph) == 2
    assert belgraph.number_of_edges() == 1

    yfg_pybel = protein('HGNC', 'PPP1R13L', identifier='18838')
    apoptosis_pybel = bioprocess('GO', name='apoptotic process',
                                 identifier='GO:0006915')
    assert yfg_pybel in belgraph
    assert apoptosis_pybel in belgraph

    _, _, e = list(belgraph.edges(data=True))[0]
    assert pc.OBJECT not in e


def test_belgraph_to_signed_graph():
    braf_no_act = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Activation(braf_no_act, mek)
    hsh = stmt.get_hash(refresh=True)

    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    pb_seg = pa.belgraph_to_signed_graph(belgraph, propagate_annotations=True)

    assert len(pb_seg.edges) == 1

    edge = (braf_dsl, map2k1_dsl, 0)
    assert edge in pb_seg.edges

    edge_dict = pb_seg.edges.get(edge)
    assert edge_dict

    assert 'stmt_hash' in edge_dict
    assert isinstance(edge_dict['stmt_hash'], int)
    assert hsh == edge_dict['stmt_hash']

    assert 'uuid' in edge_dict
    assert isinstance(edge_dict['uuid'], str)
    assert stmt.uuid == edge_dict['uuid']

    assert 'belief' in edge_dict
    assert isinstance(edge_dict['belief'], (float, int))
    assert stmt.belief == edge_dict['belief']
