from indra.statements import Agent
from indra.sources.gnbr.processor import GnbrGeneGeneProcessor


def test_standardize_agent():
    agent = GnbrGeneGeneProcessor.standardize_agent('xxx', '673')
    assert isinstance(agent, Agent)
    assert agent.name == 'BRAF'
    assert agent.db_refs.get('TEXT') == 'xxx'
    assert agent.db_refs.get('EGID') == '673'
    assert agent.db_refs.get('HGNC') == '1097'
