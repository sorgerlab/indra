from indra.assemblers import pybel_assembler as pa

def test_simple_phosphorylation():
    braf = Agent('BRAF', db_refs={'HGNC': '1097', 'UP': 'P15056'})
    mek = Agent('MAP2K1', db_refs={'HGNC': '6840', 'UP': 'Q02750'})
    stmt = Phosphorylation(braf, mek, 'S', '218')
    pba = pa.PybelAssembler([stmt])
    belgraph = pba.make_model()
    belgraph.to_bel('test.bel')
