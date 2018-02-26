import unittest
import json

from indra_db_api import api

class DbApiTestCase(unittest.TestCase):

    def setUp(self):
        api.app.testing = True
        self.app = api.app.test_client()

    def tearDown(self):
        pass

    def test_blank_response(self):
        """Test the response to an empty request."""
        resp = self.app.get('/statements/ ')
        assert resp.status_code == 400, 'Got clean response with no input.'

    def test_specific_query(self):
        """Test whether we can get a "fully" specified statement."""
        resp = self.app.get('/statements/ ?object=MAP2K1?subject=MAPK1'
                            '?action=Phosphorylation')
        assert resp.status_code == 200, \
            'Got error code %d.' % resp.status_code
        assert len(json.loads(resp.data.decode())) is not 0, \
            'Did not get any statements.'

if __name__ == '__main__':
    unittest.main()