from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import relevance_client
from indra.util import unicode_strs

small_corpus_uuid = '55c84fa4-01b4-11e5-ac0f-000c29cb28fb'

def test_get_relevant_nodes():
    nodes = relevance_client.get_relevant_nodes(small_corpus_uuid,
                                                ['MAPK1', 'MAPK3'])
    assert unicode_strs(nodes)

