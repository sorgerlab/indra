from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
from indra.databases import relevance_client
from indra.assemblers import PysbAssembler

def get_subnetwork(statements, nodes, relevance_network=None,
                   relevance_node_lim=10):
    """Return a PySB model based on a subset of given INDRA Statements.

    Statements are first filtered for nodes in the given list and other nodes
    are optionally added based on relevance in a given network. The filtered
    statements are then assembled into an executable model using INDRA's
    PySB Assembler.

    Parameters
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to extract a subnetwork from.
    nodes : list[str]
        The names of the nodes to extract the subnetwork for.
    relevance_network : Optional[str]
        The UUID of the NDEx network in which nodes relevant to the given
        nodes are found.
    relevance_node_lim : Optional[int]
        The maximal number of additional nodes to add to the subnetwork
        based on relevance.

    Returns
    -------
    model : pysb.Model
        A PySB model object assembled using INDRA's PySB Assembler from
        the INDRA Statements corresponding to the subnetwork.
    """
    if relevance_network is not None:
        relevant_nodes = _find_relevant_nodes(nodes, relevance_network,
                                              relevance_node_lim)
        all_nodes = nodes + relevant_nodes
    else:
        all_nodes = nodes
    filtered_statements = _filter_statements(statements, all_nodes)
    pa = PysbAssembler()
    pa.add_statements(filtered_statements)
    model = pa.make_model()
    return model

def _filter_statements(statements, agents):
    """Return INDRA Statements which have Agents in the given list.

    Only statements are returned in which all appearing Agents as in the
    agents list.

    Parameters
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to filter.
    agents : list[str]
        A list of agent names that need to appear in filtered statements.

    Returns
    -------
    filtered_statements : list[indra.statements.Statement]
        The list of filtered INDRA Statements.
    """
    filtered_statements = []
    for s in stmts:
        if all([a is not None for a in s.agent_list()]) and \
            all([a.name in agents for a in s.agent_list()]):
            filtered_statements.append(s)
    return filtered_statements

def _find_relevant_nodes(query_nodes, relevance_network, relevance_node_lim):
    """Return a list of nodes that are relevant for the query.

    Parameters
    ----------
    query_nodes : list[str]
        A list of node names to query for.
    relevance_network : str
        The UUID of the NDEx network to query relevance in.
    relevance_node_lim : int
        The number of top relevant nodes to return.

    Returns
    -------
    nodes : list[str]
        A list of node names that are relevant for the query.
    """
    all_nodes = relevance_client.get_relevant_nodes(relevance_network,
                                                    query_nodes)
    nodes = [n[0] for n in all_nodes[:relevance_node_lim]]
    return nodes

if __name__ == '__main__':
    genes = ['EGF', 'EGFR', 'ERBB2', 'GRB2', 'SOS1', 'HRAS', 'RAF1',
            'MAP2K1', 'MAPK1']

    with open('models/rasmachine/rem/model.pkl', 'rb') as f:
        model = pickle.load(f)
    stmts = []
    for k, v in model.items():
        stmts += v

    rasmachine_network = '50e3dff7-133e-11e6-a039-06603eb7f303'
    model = get_subnetwork(stmts, genes, rasmachine_network)

