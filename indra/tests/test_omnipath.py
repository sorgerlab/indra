import requests
from indra.sources.omnipath import OmniPathProcessor
from indra.sources.omnipath.api import op_url
from indra.statements import Agent
from indra.preassembler.grounding_mapper import GroundingMapper

BRAF_UPID = 'P15056'
JAK2_UPID = 'O60674'
CALM1_UPID = 'P0DP23'
TRPC3_UPID = 'Q13507'

BRAF_AG = Agent(None, db_refs={'UP': BRAF_UPID})
GroundingMapper.standardize_agent_name(BRAF_AG)
JAK2_AG = Agent(None, db_refs={'UP': JAK2_UPID})
GroundingMapper.standardize_agent_name(JAK2_AG)
CALM1_AG = Agent(None, db_refs={'UP': CALM1_UPID})
GroundingMapper.standardize_agent_name(CALM1_AG)
TRPC3_AG = Agent(None, db_refs={'UP': TRPC3_UPID})
GroundingMapper.standardize_agent_name(TRPC3_AG)


def test_omnipath_web_api():
    query_url = '%s/queries' % op_url
    res = requests.get(query_url)
    assert res.status_code == 200


def test_mods_from_web():
    params = {'format': 'json', 'substrates': JAK2_UPID,
              'fields': ['sources', 'references']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    assert res.status_code == 200
    assert res.text
    ptm_json = res.json()
    assert ptm_json[0]['substrate'] == JAK2_UPID, ptm_json[0]['substrate']
    op = OmniPathProcessor(ptm_json=ptm_json)
    op.process_ptm_mods()
    stmts = op.statements
    assert JAK2_AG.name in [a.name for a in stmts[0].agent_list()],\
        stmts[0].agent_list()
    assert 'omnipath' == stmts[0].evidence[0].source_api,\
        stmts[0].evidence[0].source_api


def test_ligrec_from_web():
    params = {'format': 'json', 'datasets': ['ligrecextra'],
              'fields': ['curation_effort', 'entity_type', 'references',
                         'resources', 'sources', 'type'],
              'sources': [CALM1_UPID]}
    query_url = '%s/interactions' % op_url
    res = requests.get(query_url, params)
    assert res.status_code == 200
    assert res.text
    assert 'error' not in res.text.lower()
    ligrec_json = res.json()
    assert ligrec_json[0]['source'] == CALM1_UPID
    op = OmniPathProcessor(ligrec_json=ligrec_json)
    op.process_ligrec_interactions()
    stmts = op.statements
    assert CALM1_AG.name in [a.name for a in stmts[0].agent_list()], \
        stmts[0].agent_list()
    assert 'omnipath' == stmts[0].evidence[0].source_api,\
        stmts[0].evidence[0].source_api
