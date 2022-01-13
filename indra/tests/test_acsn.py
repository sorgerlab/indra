import pandas
import requests
from indra.sources.acsn import api
from indra.sources.acsn import processor
from indra.ontology.standardize import get_standard_agent

# Agents
VEGFA = get_standard_agent('VEGFA', db_refs={'HGNC': '12680'})
MIRLET7A = get_standard_agent('MIRLET7A', db_refs={'FPLX': 'MIRLET7A'})


relations_df = pandas.read_csv(api.ACSN_RELATIONS_URL, sep='\t')
correspondence_dict = api._transform_gmt(
    requests.get(api.ACSN_CORRESPONDENCE_URL).text.split('\n'))
acsn_proc = processor.AcsnProcessor(relations_df, correspondence_dict)


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
    fplx_lookup = processor._make_famplex_lookup()
    assert 'USPL' in fplx_lookup[('CYLD', 'USPL1')]
    assert 'VEGFRR' not in fplx_lookup[('FLT1', 'FLT4', 'KDR')]


def test_get_stmt_type():
    assert processor.get_stmt_type('CATALYSIS').__name__ == 'Activation'
    assert processor.get_stmt_type('INHIBITION').__name__ == 'Inhibition'
    assert processor.get_stmt_type('HETERODIMER_ASSOCIATION').__name__ == 'Complex'
    assert processor.get_stmt_type('CATALYSIS;HETERODIMER_ASSOCIATION').__name__ == \
           'Complex'
    assert not processor.get_stmt_type('Activation')


def test_get_agent():
    assert acsn_proc.get_agent('VEGFA').db_refs == VEGFA.db_refs, VEGFA.db_refs
    assert acsn_proc.get_agent('MIRLET7A*').db_refs == MIRLET7A.db_refs, MIRLET7A.db_refs
    assert acsn_proc.get_agent('XyZ') == None

