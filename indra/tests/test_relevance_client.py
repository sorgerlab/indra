from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import relevance_client
from indra.util import unicode_strs

network_uuid = '50e3dff7-133e-11e6-a039-06603eb7f303'

def test_get_heat_kernel():
    kernel_id = relevance_client.get_heat_kernel(network_uuid)
    assert kernel_id is not None

def test_get_relevant_nodes():
    nodes = relevance_client.get_relevant_nodes(network_uuid,
                                                ['MAPK1', 'MAPK3'])
    assert unicode_strs(nodes)

