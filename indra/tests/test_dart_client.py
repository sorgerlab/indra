import json
from indra.literature.dart_client import jsonify_query_data


def test_timestamp():
    # Should ignore "after"
    assert jsonify_query_data(timestamp={'on': '2020-01-01',
                                         'after': '2020-01-02'}) == \
        json.dumps({"timestamp": {"on": "2020-01-01"}})
    assert jsonify_query_data(timestamp={'after': '2020-01-01',
                                         'before': '2020-01-05'}) == \
        json.dumps(
            {'timestamp': {'after': '2020-01-01', 'before': '2020-01-05'}})


def test_lists():
    # Check lists, ignore the lists that have non-str objects
    assert jsonify_query_data(readers=['hume', 123456],
                              versions=['123', '456']) ==\
        json.dumps({'versions': ['123', '456']})
