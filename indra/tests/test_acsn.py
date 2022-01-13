import pandas
import requests
from indra.sources.acsn import api
from indra.ontology.standardize import get_standard_agent
from indra.sources.acsn.processor import get_stmt_type, AcsnProcessor


# Agents
VEGFA = get_standard_agent('VEGFA', db_refs={'HGNC': '12680'})
MIRLET7A = get_standard_agent('MIRLET7A', db_refs={'FPLX': 'MIRLET7A'})


relations_df = pandas.read_csv(api.ACSN_RELATIONS_URL, sep='\t')
correspondence_dict = api._transform_gmt(
    requests.get(api.ACSN_CORRESPONDENCE_URL).text.split('\n'))
ap = AcsnProcessor(relations_df, correspondence_dict)


def test_acsn_web_api():
    rel_res = requests.get(api.ACSN_RELATIONS_URL)
    assert rel_res.status_code == 200

    corr_res = requests.get(api.ACSN_CORRESPONDENCE_URL)
    assert corr_res.status_code == 200


def test_transform_gmt():
    gmt_file = requests.get(api.ACSN_CORRESPONDENCE_URL).text.split('\n')
    gmt_dict = api._transform_gmt(gmt_file)
    assert 'C3' in gmt_dict['C3B*']
    assert 'ZO4*' not in gmt_dict
    assert not gmt_dict['SLC2A1'][0].endswith('\t')
    assert not gmt_dict['SLC2A1'][0].startswith('\t')


def test_famplex_lookup():
    fplx_lookup = ap.fplx_lookup
    assert 'USPL' in fplx_lookup[('CYLD', 'USPL1')]
    assert 'VEGFRR' not in fplx_lookup[('FLT1', 'FLT4', 'KDR')]


def test_get_stmt_type():
    assert get_stmt_type('CATALYSIS').__name__ == 'Activation'
    assert get_stmt_type('INHIBITION').__name__ == 'Inhibition'
    assert get_stmt_type('HETERODIMER_ASSOCIATION').__name__ == 'Complex'
    assert get_stmt_type('CATALYSIS;HETERODIMER_ASSOCIATION').__name__ == \
           'Complex'
    assert not get_stmt_type('Activation')


def test_get_agent():
    assert ap.get_agent('VEGFA').db_refs == VEGFA.db_refs, VEGFA.db_refs
    assert ap.get_agent('MIRLET7A*').db_refs == MIRLET7A.db_refs, MIRLET7A.db_refs
    assert ap.get_agent('XyZ') is None


def test_extract_statements():
    ap.extract_statements()
    stmts = ap.statements
    assert stmts[345].evidence[0].source_api == 'acsn'
    test_stmt = [stmt for stmt in stmts if(any(ag.name == 'SIVA1' for ag in stmt.agent_list()) and
                                           any(ag.name == 'TRAF2' for ag in stmt.agent_list()))]
    assert test_stmt[0] in stmts
    assert '19392652' in test_stmt[0].evidence[0].pmid
