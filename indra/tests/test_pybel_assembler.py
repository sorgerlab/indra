import networkx as nx
import pybel
import pybel.constants as pc
from indra.statements import *
from indra.databases import hgnc_client
from indra.assemblers import pybel_assembler as pa
from indra.preassembler.hierarchy_manager import hierarchies

def id(gene_name):
    return hgnc_client.get_hgnc_id(gene_name)

phostuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Ser', 218)
ubtuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ub'), 'Ser', 218)
braf_node = (pc.PROTEIN, 'HGNC', 'BRAF')
map2k1_node = (pc.PROTEIN, 'HGNC', 'MAP2K1')
tp53_node = (pc.PROTEIN, 'HGNC', 'TP53')
mdm2_node = (pc.PROTEIN, 'HGNC', 'MDM2')
egfr_node = (pc.PROTEIN, 'HGNC', 'EGFR')
egfr_phostuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Tyr', 1173)

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

def get_edge_data(g, u, v):
    return list(g.get_edge_data(u, v).values())[0]

def test_simple_modification_no_evidence():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_cat = Agent('BRAF', activity=ActivityCondition('catalytic', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt1 = Phosphorylation(braf, mek, 'S', '218')
    stmt2 = Phosphorylation(braf_kin, mek, 'S', '218')
    stmt3 = Ubiquitination(braf_cat, mek, 'S', '218')
    # Edge info for subject
    edge1 = None
    edge2 = {pc.MODIFIER: pc.ACTIVITY,
             pc.EFFECT: {
                 pc.NAME: 'kin',
                 pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}
    edge3 = {pc.MODIFIER: pc.ACTIVITY,
             pc.EFFECT: {
                 pc.NAME: 'cat',
                 pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}
    for stmt, modtuple, subj_edge in ((stmt1, phostuple, edge1),
                                     (stmt2, phostuple, edge2),
                                     (stmt3, ubtuple, edge3)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert len(belgraph.nodes()) == 3
        assert braf_node in belgraph
        map2k1_mod_node = map2k1_node + tuple([modtuple])
        assert map2k1_mod_node in belgraph
        assert belgraph.node[braf_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'BRAF'}
        assert belgraph.node[map2k1_mod_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'MAP2K1',
                pc.VARIANTS: [{
                    pc.KIND: pc.PMOD,
                    pc.IDENTIFIER: {
                        pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                        pc.NAME: modtuple[1][1]},
                    pc.PMOD_CODE: modtuple[2],
                    pc.PMOD_POSITION: modtuple[3]}]}
        assert belgraph.number_of_edges() == 2
        edge_data = get_edge_data(belgraph, braf_node, map2k1_mod_node)
        assert edge_data.get(pc.SUBJECT) == subj_edge
        assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES


def test_modification_with_mutation():
    braf = Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')],
                 db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # Adds in the base protein nodes as well as the variants (so 4 nodes)
    assert len(belgraph.nodes()) == 4
    braf_mut_node = braf_node + ((pc.HGVS, 'p.Val600Glu'),)
    assert braf_mut_node in belgraph
    assert belgraph.node[braf_mut_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'BRAF',
                pc.VARIANTS: [{
                    pc.KIND: pc.HGVS,
                    pc.IDENTIFIER: 'p.Val600Glu'}]}


def test_activation():
    braf_no_act = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt1 = Activation(braf_no_act, mek)
    stmt2 = Activation(braf_kin, mek, 'kinase')
    edge1 = {pc.RELATION: pc.DIRECTLY_INCREASES,
             pc.OBJECT: {pc.MODIFIER: pc.ACTIVITY}}
    edge2 = {pc.RELATION: pc.DIRECTLY_INCREASES,
             pc.SUBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'kin',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}},
             pc.OBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'kin',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}}
    for stmt, edge in ((stmt1, edge1), (stmt2, edge2)):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert len(belgraph.nodes()) == 2
        assert braf_node in belgraph
        assert belgraph.node[braf_node] == {
                    pc.FUNCTION: pc.PROTEIN,
                    pc.NAMESPACE: 'HGNC',
                    pc.NAME: 'BRAF'}
        assert belgraph.node[map2k1_node] == {
                    pc.FUNCTION: pc.PROTEIN,
                    pc.NAMESPACE: 'HGNC',
                    pc.NAME: 'MAP2K1'}
        assert belgraph.number_of_edges() == 1
        _, _, edge_data = belgraph.edges(data=True)[0]
        assert edge_data == edge


def test_inhibition():
    braf_kin = Agent('BRAF', activity=ActivityCondition('kinase', True),
                     db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Inhibition(braf_kin, mek, 'kinase')
    edge = {pc.RELATION: pc.DIRECTLY_DECREASES,
             pc.SUBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'kin',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}},
             pc.OBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'kin',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}}
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph.nodes()) == 2
    assert braf_node in belgraph
    assert belgraph.node[braf_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'BRAF'}
    assert belgraph.node[map2k1_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'MAP2K1'}
    assert belgraph.number_of_edges() == 1
    _, _, edge_data = belgraph.edges(data=True)[0]
    assert edge_data == edge


def test_increase_amount():
    tp53 = Agent('TP53', db_refs={'HGNC': '11998'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '6973'})

    stmt = IncreaseAmount(tp53, mdm2)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph.nodes()) == 2
    assert mdm2_node in belgraph
    assert tp53_node in belgraph
    assert belgraph.node[tp53_node] == {pc.FUNCTION: pc.PROTEIN,
                                        pc.NAMESPACE: 'HGNC',
                                        pc.NAME: 'TP53'}
    assert belgraph.node[mdm2_node] == {pc.FUNCTION: pc.PROTEIN,
                                        pc.NAMESPACE: 'HGNC',
                                        pc.NAME: 'MDM2'}
    assert belgraph.number_of_edges() == 1
    _, _, edge_data = belgraph.edges(data=True)[0]
    assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES


def test_increase_amount_tscript():
    tp53 = Agent('TP53', activity=ActivityCondition('transcription', True),
                 db_refs={'HGNC': '11998'})
    mdm2 = Agent('MDM2', db_refs={'HGNC': '6973'})

    stmt = IncreaseAmount(tp53, mdm2)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph.nodes()) == 2
    assert mdm2_node in belgraph
    assert tp53_node in belgraph
    assert belgraph.node[tp53_node] == {pc.FUNCTION: pc.PROTEIN,
                                        pc.NAMESPACE: 'HGNC',
                                        pc.NAME: 'TP53'}
    assert belgraph.node[mdm2_node] == {pc.FUNCTION: pc.PROTEIN,
                                        pc.NAMESPACE: 'HGNC',
                                        pc.NAME: 'MDM2'}
    assert belgraph.number_of_edges() == 1
    _, _, edge_data = belgraph.edges(data=True)[0]
    assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES
    assert edge_data[pc.SUBJECT] == {
            pc.MODIFIER: pc.ACTIVITY,
            pc.EFFECT: {pc.NAME: 'tscript',
                        pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}


def test_gef():
    gef = Agent('SOS1', mods=[ModCondition('phosphorylation')],
                db_refs={'HGNC':'11187'})
    ras = Agent('KRAS', db_refs={'HGNC':'6407'})
    stmt = Gef(gef, ras)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 3
    assert belgraph.number_of_edges() == 2
    gef_node = (pc.PROTEIN, 'HGNC', 'SOS1',
                (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph')))
    kras_node = (pc.PROTEIN, 'HGNC', 'KRAS')
    assert gef_node in belgraph
    assert kras_node in belgraph
    edge_data = get_edge_data(belgraph, gef_node, kras_node)
    edge = {pc.RELATION: pc.DIRECTLY_INCREASES,
             pc.SUBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'gef',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}},
             pc.OBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'gtp',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}}
    assert edge_data == edge


def test_gap():
    gap = Agent('RASA1', mods=[ModCondition('phosphorylation')],
                db_refs={'HGNC':'9871'})
    ras = Agent('KRAS', db_refs={'HGNC':'6407'})
    stmt = Gap(gap, ras)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 3
    assert belgraph.number_of_edges() == 2
    gap_node = (pc.PROTEIN, 'HGNC', 'RASA1',
                 (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph')))
    ras_node = (pc.PROTEIN, 'HGNC', 'KRAS')
    assert gap_node in belgraph
    assert ras_node in belgraph
    edge_data = get_edge_data(belgraph, gap_node, ras_node)
    edge = {pc.RELATION: pc.DIRECTLY_DECREASES,
             pc.SUBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'gap',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}},
             pc.OBJECT: {
                 pc.MODIFIER: pc.ACTIVITY,
                 pc.EFFECT: {
                     pc.NAME: 'gtp',
                     pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}}
    assert edge_data == edge


def test_active_form():
    ras = Agent('KRAS', mutations=[MutCondition('12', 'G', 'V')],
                db_refs={'HGNC':'6407'})
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
    for stmt in (stmt1, stmt2, stmt3):
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        assert len(belgraph) == 2


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
    complex_tuple = [n for n in belgraph.nodes() if n[0] == pc.COMPLEX][0]
    # Don't bother to check the tuple, which is now generated by
    # PyBEL directly, but check the node data
    assert belgraph.node[complex_tuple] == {
        pc.FUNCTION: pc.COMPLEX,
        pc.MEMBERS: [
            {pc.FUNCTION: pc.PROTEIN,
             pc.NAMESPACE: 'HGNC',
             pc.NAME: 'EGFR'},
            {pc.FUNCTION: pc.PROTEIN,
             pc.NAMESPACE: 'HGNC',
             pc.NAME: 'GRB2'},
            {pc.FUNCTION: pc.PROTEIN,
             pc.NAMESPACE: 'HGNC',
             pc.NAME: 'SOS1'},
        ]}


def test_rxn_no_controller():
    glu = Agent('D-GLUCOSE', db_refs = {'CHEBI': 'CHEBI:17634'})
    g6p = Agent('GLUCOSE-6-PHOSPHATE', db_refs={'CHEBI': 'CHEBI:4170'})
    stmt = Conversion(None, [glu], [g6p])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # The graph should contain the node for the reaction as well as nodes
    # for all of the members
    assert len(belgraph) == 3
    rxn_tuple = [n for n in belgraph.nodes() if n[0] == pc.REACTION][0]
    assert belgraph.node[rxn_tuple] == {
        pc.FUNCTION: pc.REACTION,
        pc.REACTANTS: [
            {
                pc.FUNCTION: pc.ABUNDANCE,
                pc.NAMESPACE: 'CHEBI',
                pc.NAME: '17634',
            }
        ],
        pc.PRODUCTS: [
            {
                pc.FUNCTION: pc.ABUNDANCE,
                pc.NAMESPACE: 'CHEBI',
                pc.NAME: '4170'
            }
        ]
    }


def test_rxn_with_controller():
    hk1 = Agent('HK1', db_refs={'HGNC': id('HK1')})
    glu = Agent('D-GLUCOSE', db_refs={'CHEBI': 'CHEBI:17634'})
    g6p = Agent('GLUCOSE-6-PHOSPHATE', db_refs={'CHEBI': 'CHEBI:4170'})
    stmt = Conversion(hk1, [glu], [g6p])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    # The graph should contain the node for the reaction as well as nodes
    # for all of the members
    assert len(belgraph) == 4
    rxn_tuple = [n for n in belgraph.nodes() if n[0] == pc.REACTION][0]
    # The reaction data should be the same as before
    assert belgraph.node[rxn_tuple] == {
        pc.FUNCTION: pc.REACTION,
        pc.REACTANTS: [
            {
                pc.FUNCTION: pc.ABUNDANCE,
                pc.NAMESPACE: 'CHEBI',
                pc.NAME: '17634',
            }
        ],
        pc.PRODUCTS: [
            {
                pc.FUNCTION: pc.ABUNDANCE,
                pc.NAMESPACE: 'CHEBI',
                pc.NAME: '4170'
            }
        ]
    }


def test_autophosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    stmt = Autophosphorylation(egfr, 'Y', '1173')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 2
    assert egfr_node in belgraph.nodes()
    egfr_phos_node = egfr_node + tuple([egfr_phostuple])
    assert egfr_phos_node in belgraph.nodes()
    assert belgraph.node[egfr_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'EGFR'}
    assert belgraph.node[egfr_phos_node] == {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'EGFR',
                pc.VARIANTS: [{
                    pc.KIND: pc.PMOD,
                    pc.IDENTIFIER: {
                        pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                        pc.NAME: egfr_phostuple[1][1]},
                    pc.PMOD_CODE: egfr_phostuple[2],
                    pc.PMOD_POSITION: egfr_phostuple[3]}]}
    assert belgraph.number_of_edges() == 2
    # There will be two edges between these nodes
    edge_dicts = list(belgraph.get_edge_data(egfr_node,
                                             egfr_phos_node).values())
    assert {pc.RELATION: pc.DIRECTLY_INCREASES} in edge_dicts

    # Test an autophosphorylation with a bound condition
    tab1 = Agent('TAB1', db_refs={'HGNC': id('TAB1')})
    p38_tab1 = Agent('MAPK14', bound_conditions=[BoundCondition(tab1)],
                     db_refs={'HGNC': id('MAPK14')})
    stmt = Autophosphorylation(p38_tab1, 'Y', '100')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 4
    assert belgraph.number_of_edges() == 4


def test_bound_condition():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    grb2 = Agent('GRB2', db_refs={'HGNC': id('GRB2')})
    ras = Agent('KRAS', db_refs={'HGNC':'6407'})
    sos1_bound = Agent('SOS1', mods=[ModCondition('phosphorylation')],
                  bound_conditions=[BoundCondition(egfr), BoundCondition(grb2)],
                  db_refs={'HGNC': id('SOS1')})
    stmt = Gef(sos1_bound, ras)
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 6
    assert belgraph.number_of_edges() == 5
    complex_tuple = [n for n in belgraph.nodes() if n[0] == pc.COMPLEX][0]
    # Don't bother to check the tuple, which is now generated by
    # PyBEL directly, but check the node data
    assert belgraph.node[complex_tuple][pc.FUNCTION] == pc.COMPLEX
    complex_members = belgraph.node[complex_tuple][pc.MEMBERS]
    members_list = [{
                 pc.FUNCTION: pc.PROTEIN,
                 pc.NAMESPACE: 'HGNC',
                 pc.NAME: 'EGFR'
             }, {
                 pc.FUNCTION: pc.PROTEIN,
                 pc.NAMESPACE: 'HGNC',
                 pc.NAME: 'GRB2'
             }, {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'SOS1',
                pc.VARIANTS: [{
                    pc.KIND: pc.PMOD,
                    pc.IDENTIFIER: {
                        pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                        pc.NAME: 'Ph'}}]}]
    for member in members_list:
        assert member in complex_members
    kras_node = (pc.PROTEIN, 'HGNC', 'KRAS')
    assert kras_node in belgraph.nodes()
    assert (complex_tuple, kras_node) in belgraph.edges()
    edge_data = (complex_tuple, kras_node,
                {pc.RELATION: pc.DIRECTLY_INCREASES,
                 pc.OBJECT: {
                     pc.MODIFIER: pc.ACTIVITY,
                     pc.EFFECT: {
                         pc.NAME: 'gtp',
                         pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}})
    assert edge_data in belgraph.edges(data=True)


def test_transphosphorylation():
    egfr = Agent('EGFR', db_refs={'HGNC': id('EGFR')})
    egfr_dimer = Agent('EGFR', bound_conditions=[BoundCondition(egfr)],
                       db_refs={'HGNC': id('EGFR')})
    stmt = Transphosphorylation(egfr_dimer, 'Y', '1173')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph) == 3
    assert belgraph.number_of_edges() == 3
    egfr_dimer_node = (pc.COMPLEX, (pc.PROTEIN, 'HGNC', 'EGFR'),
                                   (pc.PROTEIN, 'HGNC', 'EGFR'))
    egfr_phos_node = (pc.PROTEIN, 'HGNC', 'EGFR',
                      (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Tyr', 1173))
    edge_data = get_edge_data(belgraph, egfr_dimer_node, egfr_phos_node)
    assert edge_data == {pc.RELATION: pc.DIRECTLY_INCREASES}

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
    assert len(belgraph) == 5
    members_list = [{
                 pc.FUNCTION: pc.PROTEIN,
                 pc.NAMESPACE: 'HGNC',
                 pc.NAME: 'EGFR'
             }, {
                 pc.FUNCTION: pc.PROTEIN,
                 pc.NAMESPACE: 'HGNC',
                 pc.NAME: 'GRB2'
             }, {
                pc.FUNCTION: pc.PROTEIN,
                pc.NAMESPACE: 'HGNC',
                pc.NAME: 'SOS1',
                pc.VARIANTS: [{
                    pc.KIND: pc.PMOD,
                    pc.IDENTIFIER: {
                        pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                        pc.NAME: 'Ph'},
                    pc.PMOD_CODE: 'Tyr',
                    pc.PMOD_POSITION: 100}]}]
    # Get the complex node
    cplx_node_list = [n for n in belgraph.nodes(data=True)
                 if n[0][0] == pc.COMPLEX]
    assert len(cplx_node_list) == 1
    cplx_node = cplx_node_list[0]
    cplx_node_data = cplx_node[1]
    # Check the constituents of the complex node
    for member in members_list:
        assert member in cplx_node_data[pc.MEMBERS]


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
    cplx_node_list = [n for n in belgraph.nodes() if n[0] == pc.COMPLEX]
    assert len(cplx_node_list) == 2

