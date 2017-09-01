from indra.sources.ndex_cx import process_cx_file
from indra.sources.ndex_cx.processor import NdexCxProcessor
from indra.databases import hgnc_client
from indra.statements import Agent

def test_process_cx_file():
    ncp = process_cx_file('merged_BRCA1_formatted.cx')
    assert isinstance(ncp, NdexCxProcessor)

def test_initialize_node_agents():
    ncp = process_cx_file('merged_BRCA1_formatted.cx')
    assert isinstance(ncp._node_map, dict)
    for ndex_id, agent in ncp._node_map.items():
        assert isinstance(ndex_id, int)
        assert isinstance(agent, Agent)
        #assert agent.db_refs.get('HGNC')

if __name__ == '__main__':
    #test_process_cx_file()
    #ncp = process_cx_file('merged_BRCA1_formatted.cx')
    test_initialize_node_agents()
