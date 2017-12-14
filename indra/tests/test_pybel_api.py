import pybel
from pybel.dsl import *
import pybel.constants as pc
from pybel.examples import egf_graph, sialic_acid_graph
from indra.statements import *
from indra.sources import pybel as pb
from indra.databases import hgnc_client


def test_process_pybel():
    pbp = pb.process_pybel_graph(egf_graph)
    assert pbp.statements


def test_increase_amount():
    mek = protein(name='MAP2K1', namespace='HGNC')
    erk = protein(name='MAPK1', namespace='HGNC')
    g = pybel.BELGraph()
    g.add_qualified_edge(mek, erk, relation=pc.INCREASES,
                         evidence="Some evidence.", citation='123456')
    pbp = pb.process_pybel_graph(g)
    assert pbp.statements
    assert len(pbp.statements) == 1
    assert isinstance(pbp.statements[0], IncreaseAmount)

def test_get_agent():
    mek = protein(name='MAP2K1', namespace='HGNC')
    agent = pb._get_agent(mek)
    hgnc_id = hgnc_client.get_hgnc_id('MAP2K1')
    up_id = hgnc_client.get_uniprot_id(hgnc_id)
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == hgnc_id
    assert agent.db_refs.get('UP') == up_id

    # Now create an agent with an identifier
    mek = protein(name='Foo', namespace='HGNC', identifier='6840')
    agent = pb._get_agent(mek)
    assert isinstance(agent, Agent)
    assert agent.name == 'MAP2K1'
    assert agent.db_refs.get('HGNC') == hgnc_id
    assert agent.db_refs.get('UP') == up_id


if __name__ == '__main__':
    test_get_agent()
