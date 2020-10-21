import json
import requests
from nose.plugins.attrib import attr
from indra.config import get_config
from indra.literature import dart_client


def test_timestamp():
    # Should ignore "after"
    assert dart_client._jsonify_query_data(
        timestamp={'on': '2020-01-01',
                   'after': '2020-01-02'}) == \
        json.dumps({"timestamp": {"on": "2020-01-01"}})
    assert dart_client._jsonify_query_data(
        timestamp={'after': '2020-01-01',
                   'before': '2020-01-05'}) == \
        json.dumps(
            {'timestamp': {'after': '2020-01-01', 'before': '2020-01-05'}})


def test_lists():
    # Check lists, ignore the lists that have non-str objects
    assert dart_client._jsonify_query_data(readers=['hume', 123456],
                                           versions=['123', '456']) == \
        json.dumps({'versions': ['123', '456']})


def test_prioritize():
    records = \
        [{'identity': 'cwms',
          'version': '20200819',
          'document_id': '5cd8e79aef0d08cd9ccb758f5157ec22',
          'storage_key': '398387b5-7c71-41a0-9d2a-eee6d4e5ed59.ekb'},
         {'identity': 'cwms',
          'version': '20200819',
          'document_id': '8ef92f70c7ff31b75883aaac49fb71e5',
          'storage_key': '46250fd9-5ac8-4b8a-931b-62bf77bcd154.ekb'},
         {'identity': 'cwms',
          'version': '20200819',
          'document_id': 'c5fe18952cffe87fe77c271a1ab92130',
          'storage_key': '25aecbdf-a632-41ff-80a2-883f7f47a830.ekb'},
         {'identity': 'cwms',
          'version': '20200820',
          'document_id': '1f1f612a88418c6c2ca3fe88aeae459e',
          'storage_key': '81eb10b5-2c71-4a6f-8553-f3e3ae3c1147.ekb'},
         {'identity': 'cwms',
          'version': '20200820',
          'document_id': '51beb529d3b14fe93e5febacc7d39b91',
          'storage_key': '9d6029d0-fc20-41f0-a736-ad03c1421cad.ekb'},
         {'identity': 'cwms',
          'version': '20200821',
          'document_id': '1f1f612a88418c6c2ca3fe88aeae459e',
          'storage_key': '6f6fdf90-ff4f-45db-93b7-d8414e842cb3.ekb'},
         {'identity': 'cwms',
          'version': '20200821',
          'document_id': '51beb529d3b14fe93e5febacc7d39b91',
          'storage_key': 'c7bbaead-1bac-43cb-aa5a-d934baad145a.ekb'},
         {'identity': 'cwms',
          'version': '20200821',
          'document_id': '5cd8e79aef0d08cd9ccb758f5157ec22',
          'storage_key': '322cd716-b13a-4477-b005-9cb7d71df031.ekb'},
         {'identity': 'cwms',
          'version': '20200821',
          'document_id': '8ef92f70c7ff31b75883aaac49fb71e5',
          'storage_key': '71c956f9-e6c3-4648-a007-3e30411be206.ekb'},
         {'identity': 'cwms',
          'version': '20200821',
          'document_id': 'c5fe18952cffe87fe77c271a1ab92130',
          'storage_key': 'c0786283-6f5d-436d-a8a2-d83a6e98de60.ekb'}
         ]

    prioritized_records = dart_client.prioritize_records(
        records,
        priorities={'cwms': ['20200821', '20200820', '20200819']})
    assert len(prioritized_records) == 5, prioritized_records
    assert all([rec['version'] == '20200821' for rec in prioritized_records])

    prioritized_records = dart_client.prioritize_records(
        records,
        priorities={'cwms': ['20200819', '20200820', '20200821']})
    assert len(prioritized_records) == 5, prioritized_records
    assert all([rec['version'] != '20200821' for rec in prioritized_records])


@attr('nonpublic', 'notravis')
def test_api():
    health_ep = dart_client.dart_base_url + '/health'
    dart_uname = get_config('DART_WM_USERNAME', failure_ok=False)
    dart_pwd = get_config('DART_WM_PASSWORD', failure_ok=False)
    res = requests.get(health_ep, auth=(dart_uname, dart_pwd))
    assert res.status_code == 200
