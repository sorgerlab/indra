from indra.assemblers import pybel_assembler as pa
from indra.statements import *
import pybel
import networkx as nx
import pybel.constants as pc

phostuple = (pc.PMOD, (pc.BEL_DEFAULT_NAMESPACE, 'Ph'), 'Ser', 218)

def draw(g, filename):
    ag = nx.nx_agraph.to_agraph(g)
    ag.draw(filename, prog='dot')

def test_simple_phosphorylation_no_evidence():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    #pybel.to_bel(belgraph)
    assert len(belgraph.nodes()) == 2
    assert ('Protein', 'HGNC', 'BRAF') in belgraph
    assert ('Protein', 'HGNC', 'MAP2K1', phostuple) in belgraph
    assert len(belgraph.edges()) == 1
    _, _, edge_data = belgraph.edges(data=True)[0]
    assert edge_data[pc.SUBJECT] == {
                pc.MODIFIER: pc.ACTIVITY,
                pc.EFFECT: {
                    pc.NAME: 'kin', pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}}
    assert edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES

    draw(belgraph, 'belgraph.pdf')

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

