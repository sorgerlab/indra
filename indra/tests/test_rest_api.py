import requests
from nose.plugins.attrib import attr

@attr('webservice')
def test_rest_api_responsive():
    stmt_str = '{"statements": [{"sbo": "http://identifiers.org/sbo/SBO:0000526", "type": "Complex", "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f", "members": [{"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"}, {"db_refs": {"TEXT": "ERK", "NCIT": "C26360", "FPLX": "ERK"}, "name": "ERK"}], "evidence": [{"text": "MEK binds ERK", "source_api": "trips"}]}]}'
    url = 'http://api.indra.bio:8000/' + \
        'assemblers/cyjs'
    res = requests.post(url, stmt_str)
    assert res.status_code == 200
