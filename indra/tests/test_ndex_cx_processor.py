from indra.sources.ndex_cx import process_cx_file, process_ndex_network
from indra.sources.ndex_cx.processor import NdexCxProcessor
from indra.databases import hgnc_client
from indra.statements import Agent, Statement
from requests.exceptions import HTTPError
from nose.tools import raises


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


def test_get_pmids():
    pmids = ncp_file.get_pmids()
    for pmid in pmids:
        # Make sure all of the entries are valid integers
        assert isinstance(pmid, str)
        int(pmid)


def test_get_statements():
    stmts = ncp_file.get_statements()
    for stmt in stmts:
        assert isinstance(stmt, Statement)
        for ag in stmt.agent_list():
            assert isinstance(ag, Agent)
        for ev in stmt.evidence:
            assert ev.source_api == 'ndex'


def test_get_cx_from_ndex():
    # This network is public
    ncp = process_ndex_network('171b8e16-8cf4-11e7-a10d-0ac135e8bacf')


@raises(HTTPError)
def test_get_cx_from_ndex_unauth():
    # This network should error because unauthorized without username/pwd
    ncp = process_ndex_network('df1fea48-8cfb-11e7-a10d-0ac135e8bacf')
