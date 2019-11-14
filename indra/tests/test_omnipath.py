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

def test_query_ptms():
    stmts = op.get_ptms(['Q13873'])
    assert len(stmts) == 1
    assert isinstance(stmts[0], Phosphorylation)
    assert stmts[0].enz.name == 'CSNK2A1'
    assert stmts[0].sub.name == 'BMPR2'
    assert stmts[0].residue == 'S'
    assert stmts[0].position == '757'
