import json
import requests
from indra.config import get_config
from indra.literature.dart_client import _jsonify_query_data, dart_base_url


def test_timestamp():
    # Should ignore "after"
    assert _jsonify_query_data(timestamp={'on': '2020-01-01',
                                          'after': '2020-01-02'}) == \
        json.dumps({"timestamp": {"on": "2020-01-01"}})
    assert _jsonify_query_data(timestamp={'after': '2020-01-01',
                                          'before': '2020-01-05'}) == \
        json.dumps(
            {'timestamp': {'after': '2020-01-01', 'before': '2020-01-05'}})


def test_lists():
    # Check lists, ignore the lists that have non-str objects
    assert _jsonify_query_data(readers=['hume', 123456],
                               versions=['123', '456']) ==\
        json.dumps({'versions': ['123', '456']})


def test_api():
    health_ep = dart_base_url + '/health'
    dart_uname = get_config('DART_WM_USERNAME', failure_ok=False)
    dart_pwd = get_config('DART_WM_PASSWORD', failure_ok=False)
    res = requests.get(health_ep, auth=(dart_uname, dart_pwd))
    assert res.status_code == 200
