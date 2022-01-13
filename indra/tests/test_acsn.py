import requests
from indra.sources.acsn import api
from indra.sources.acsn import processor


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
