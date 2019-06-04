import json
import requests
from datetime import datetime

from nose.plugins.attrib import attr

BASE = 'http://api.indra.bio:8000/'


def _call_api(method, route, *args, **kwargs):
    route = route.lstrip('/')
    req_meth = getattr(requests, method)
    start = datetime.now()
    print("Submitting request to '%s' at %s." % ('/' + route, start))
    print("\targs:", args)
    print("\tkwargs:", kwargs)
    res = req_meth(BASE + route, *args, **kwargs)
    end = datetime.now()
    print("Got result with %s at %s after %s seconds."
          % (res.status_code, end, (end-start).total_seconds()))
    assert res.status_code == 200, res.status_code
    return res


@attr('webservice')
def test_responsive():
    res = _call_api('get', '')
    assert res.content.startswith(b'This is the INDRA REST API.')


@attr('webservice')
def test_assemblers_cyjs():
    stmt_json = {
        "statements": [{
            "sbo": "http://identifiers.org/sbo/SBO:0000526",
            "type": "Complex",
            "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
            "members": [
                {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
                {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"}
            ],
            "evidence": [{"text": "MEK binds ERK", "source_api": "trips"}]
        }]
    }
    stmt_str = json.dumps(stmt_json)
    res = _call_api('post', 'assemblers/cyjs', stmt_str)
    res_json = res.json()
    assert len(res_json['edges']) == 1, len(res_json['edges'])
    assert len(res_json['nodes']) == 2, len(res_json['nodes'])
    return
