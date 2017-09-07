from indra.assemblers import pybel_assembler as pa
from indra.statements import *
import pybel
import networkx as nx
import pybel.constants as pc

phostuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Ser', 218)
ubtuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ub'), 'Ser', 218)
braf_node = (pc.PROTEIN, 'HGNC', 'BRAF')
map2k1_node = (pc.PROTEIN, 'HGNC', 'MAP2K1')
tp53_node = (pc.PROTEIN, 'HGNC', 'TP53')
mdm2_node = (pc.PROTEIN, 'HGNC', 'MDM2')

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

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
        assert len(belgraph.nodes()) == 2
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
        assert belgraph.number_of_edges() == 1
        _, _, edge_data = belgraph.edges(data=True)[0]
        assert edge_data.get(pc.SUBJECT) == subj_edge
        assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES


def test_modification_with_mutation():
    braf = Agent('BRAF', mutations=[MutCondition('600', 'V', 'E')],
                 db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    assert len(belgraph.nodes()) == 2
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

