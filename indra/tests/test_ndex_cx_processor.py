from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.sources.ndex_cx import process_cx_file, process_ndex_network
from indra.sources.ndex_cx.processor import NdexCxProcessor
from indra.statements import Agent, Statement
from requests.exceptions import HTTPError
import pytest


path_this = os.path.dirname(os.path.abspath(__file__))
ncp_file = \
    process_cx_file(os.path.join(path_this, 'merged_BRCA1_formatted.cx'))


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


def test_get_agents_no_grounding():
    ncp = \
        process_cx_file(os.path.join(path_this, 'merged_BRCA1_formatted.cx'),
                        require_grounding=False)
    node_names = list(ncp._node_names.values())
    names_from_agents = [ag.name for ag in ncp._node_agents.values()]
    texts_from_agents = [ag.db_refs['TEXT'] for ag in ncp._node_agents.values()]
    assert set(node_names) == set(names_from_agents)
    assert set(node_names) == set(texts_from_agents)


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


@pytest.mark.webservice
def test_get_cx_from_ndex():
    # Ras Machine network
    ncp = process_ndex_network('fc56fe8d-1b60-11e8-b939-0ac135e8bacf')


@pytest.mark.webservice
def test_get_cx_from_ndex_unauth():
    # This network should error because unauthorized without username/pwd
    with pytest.raises(HTTPError):
        ncp = process_ndex_network('df1fea48-8cfb-11e7-a10d-0ac135e8bacf')
