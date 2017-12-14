import pybel
from pybel.examples import egf_graph, sialic_acid_graph
from indra.sources import pybel as pb

def test_process_pybel():
    pbp = pb.process_pybel_graph(egf_graph)
    assert pbp.statements

