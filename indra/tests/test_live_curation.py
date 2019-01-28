import json
import unittest
from indra.tools.live_curation import app


class LiveCurationTestCase(unittest.TestCase):

    def setUp(self):
        app.testing = True
        self.app = app.test_client()

    # Tests ==================
    def test_alive(self):
        resp = self.app.post('update_beliefs',
                             data=json.dumps({'content': 'nothing'}))
        assert resp.status_code != 200, 'Bogus post worked!'
