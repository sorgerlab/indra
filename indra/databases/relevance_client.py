from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from indra.databases import ndex_client
# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('relevance')

ndex_relevance = 'http://general.bigmech.ndexbio.org:8080'

def get_heat_kernel(network_id):
    """Return the identifier of a heat kernel calculated for a given network.

    Parameters
    ----------
    network_id : str
        The UUID of the network in NDEx.

    Returns
    -------
    kernel_id : str
        The identifier of the heat kernel calculated for the given network.
    """
    url = ndex_relevance + '/%s/generate_ndex_heat_kernel' % network_id
    res = ndex_client.send_request(url, {}, is_json=True, use_get=True)
    kernel_id = res.get('kernel_id')
    if kernel_id is None:
        logger.error('Could not get heat kernel for network.')
        return None
    return kernel_id

def get_relevant_nodes(network_id, query_nodes):
    """Return a set of network nodes relevant to a given query set.

    A heat diffusion algorithm is used on a pre-computed heat kernel for the
    given network which starts from the given query nodes. The nodes
    in the network are ranked according to heat score which is a measure
    of relevance with respect to the query nodes.

    Parameters
    ----------
    network_id : str
        The UUID of the network in NDEx.
    query_nodes : list[str]
        A list of node names with respect to which relevance is queried.

    Returns
    -------
    ranked_entities : list[(str, float)]
        A list containing pairs of node names and their relevance scores.
    """
    url = ndex_relevance + '/rank_entities'
    kernel_id = get_heat_kernel(network_id)
    if isinstance(query_nodes, basestring):
        query_nodes = [query_nodes]
    params = {'identifier_set': query_nodes,
              'kernel_id': kernel_id}
    res = ndex_client.send_request(url, params, is_json=True)
    ranked_entities = res.get('ranked_entities')
    if ranked_entities is None:
        logger.error('Could not get ranked entities.')
        return None
    return ranked_entities
