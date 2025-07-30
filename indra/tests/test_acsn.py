import pandas
import requests
from unittest import skip
from functools import lru_cache
from indra.sources.acsn import api
from indra.ontology.standardize import get_standard_agent
from indra.sources.acsn.processor import get_stmt_type, AcsnProcessor


@lru_cache(1)
def _get_relations_df():
    relations_df = pandas.read_csv(api.ACSN_RELATIONS_URL, sep='\t')
    return relations_df


@lru_cache(1)
def _get_gmt_file():
    gmt_file = requests.get(api.ACSN_CORRESPONDENCE_URL).text.split('\n')
    return gmt_file


@lru_cache(1)
def _get_acsn_processor():
    gmt_file = _get_gmt_file()
    relations_df = _get_relations_df()
    correspondence_dict = api._transform_gmt(gmt_file)
    ap = AcsnProcessor(relations_df, correspondence_dict)
    return ap


@skip("https://acsn.curie.fr/ is down")
def test_transform_gmt():
    gmt_file = _get_gmt_file()
    gmt_dict = api._transform_gmt(gmt_file)
    assert 'C3' in gmt_dict['C3B*']
    assert 'ZO4*' not in gmt_dict
    assert not gmt_dict['SLC2A1'][0].endswith('\t')
    assert not gmt_dict['SLC2A1'][0].startswith('\t')


@skip("https://acsn.curie.fr/ is down")
def test_famplex_lookup():
    ap = _get_acsn_processor()
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


@skip("https://acsn.curie.fr/ is down")
def test_get_agent():
    ap = _get_acsn_processor()
    # Agents
    VEGFA = get_standard_agent('VEGFA', db_refs={'HGNC': '12680'})
    MIRLET7A = get_standard_agent('MIRLET7A', db_refs={'FPLX': 'MIRLET7A'})

    assert ap.get_agent('VEGFA').db_refs == VEGFA.db_refs, VEGFA.db_refs
    assert ap.get_agent('MIRLET7A*').db_refs == \
           MIRLET7A.db_refs, MIRLET7A.db_refs
    assert ap.get_agent('XyZ') is None


@skip("https://acsn.curie.fr/ is down")
def test_extract_statements():
    ap = _get_acsn_processor()
    ap.extract_statements()
    stmts = ap.statements
    assert stmts[345].evidence[0].source_api == 'acsn'
    test_stmt = [stmt for stmt in stmts if(any(ag.name == 'SIVA1'
                                               for ag in stmt.agent_list()) and
                                           any(ag.name == 'TRAF2'
                                               for ag in stmt.agent_list()))]
    assert test_stmt[0] in stmts
    assert '19392652' in test_stmt[0].evidence[0].pmid
