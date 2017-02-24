from indra.assemblers import pybel_assembler as pa
from indra.statements import *
import pybel

def test_simple_phosphorylation():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    pybel.to_bel(belgraph)
    assert len(belgraph.nodes()) == 2
    assert ('Protein', 'HGNC', 'BRAF') in belgraph.nodes()
    assert ('Protein', 'HGNC', 'MAP2K1') in belgraph.nodes()
    assert len(belgraph.edges()) == 1

if __name__ == '__main__':
     test_simple_phosphorylation()

