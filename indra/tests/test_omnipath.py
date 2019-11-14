import requests
from indra.sources.omnipath import OmniPathModificationProcessor,\
    OmniPathLiganReceptorProcessor
from indra.sources.omnipath.api import op_url
from indra.statements import Agent, Phosphorylation
from indra.preassembler.grounding_mapper import GroundingMapper

BRAF_UPID = 'P15056'
JAK2_UPID = 'O60674'
BRAF_AG = Agent(None, db_refs={'UP': BRAF_UPID})
GroundingMapper.standardize_agent_name(BRAF_AG)
JAK2_AG = Agent(None, db_refs={'UP': JAK2_UPID})
GroundingMapper.standardize_agent_name(JAK2_AG)


def test_omnipath_web_api():
    query_url = '%s/queries'
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
    stmts = OmniPathModificationProcessor(ptm_json).statements
    assert JAK2_AG.name in [a.name for a in stmts[0].agent_list()],\
        stmts[0].agent_list()
    assert 'omnipath' == stmts[0].evidence[0].source_api,\
        stmts[0].evidence[0].source_api
