from indra.sources.ndex_cx import process_cx_file
from indra.sources.ndex_cx.processor import NdexCxProcessor
from indra.databases import hgnc_client
from indra.statements import Agent, Statement

ncp_file = process_cx_file('merged_BRCA1_formatted.cx')

def test_process_cx_file():
    assert isinstance(ncp_file, NdexCxProcessor)

def test_initialize_node_agents():
    assert isinstance(ncp_file._node_agents, dict)
    for ndex_id, agent in ncp_file._node_agents.items():
        assert isinstance(ndex_id, int)
        assert isinstance(agent, Agent)
        assert agent.db_refs.get('HGNC')
        assert agent.db_refs.get('UP')

def test_get_agents():
    nodes = ncp_file.get_agents()
    assert nodes == list(ncp_file._node_agents.values())

def test_get_node_names():
    nodes = ncp_file.get_node_names()
    assert nodes == list(ncp_file._node_names.values())

def test_get_statements():
    stmts = ncp_file.get_statements()
    for stmt in stmts:
        assert isinstance(stmt, Statement)
        for ag in stmt.agent_list():
            assert isinstance(ag, Agent)

if __name__ == '__main__':
    #test_process_cx_file()
    #ncp = process_cx_file('merged_BRCA1_formatted.cx')
    #import ipdb; ipdb.set_trace()
    #test_initialize_node_agents()
    #stmts = test_get_statements()
    pass
