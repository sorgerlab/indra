from indra.assemblers import pybel_assembler as pa
from indra.statements import *
import pybel
import networkx as nx
import pybel.constants as pc

phostuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Ser', 218)
ubtuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ub'), 'Lys', 218)
braf_node = (pc.PROTEIN, 'HGNC', 'BRAF')
map2k1_node = (pc.PROTEIN, 'HGNC', 'MAP2K1')

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

def test_simple_phosphorylation_no_evidence():
    for modclass, modtuple in ((Phosphorylation, phostuple),
                               (Ubiquitination, ubtuple)):
        braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
        mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
        stmt = modclass(braf, mek, 'S', '218')
        pba = pa.PybelAssembler([stmt])
        belgraph = pba.make_model()
        import ipdb; ipdb.set_trace()
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
        assert edge_data[pc.SUBJECT] == {
                pc.MODIFIER: pc.ACTIVITY,
                pc.EFFECT: {
                    pc.NAME: 'kin', pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}
        assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES

def test_simple_phosphorylation_with_evidence():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    ev1 = Evidence('reach', '1234', '12345678', 'Some example text.')
    ev2 = Evidence('reach', '2345', '23456789', 'Some more example text.')
    stmt = Phosphorylation(braf, mek, 'S', '218', evidence=[ev1, ev2])
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    pybel.to_bel(belgraph)
    assert len(belgraph.nodes()) == 2
    assert ('Protein', 'HGNC', 'BRAF') in belgraph.nodes()
    assert ('Protein', 'HGNC', 'MAP2K1') in belgraph.nodes()
    assert len(belgraph.edges()) == 2

if __name__ == '__main__':
     test_simple_phosphorylation_no_evidence()

