import json
import logging
from os import path

import numpy as np
import pandas as pd
import networkx as nx
from decimal import Decimal

import indra
from indra.belief import SimpleScorer
from indra.statements import Evidence
from indra.statements import Statement

logger = logging.getLogger(__name__)
simple_scorer = SimpleScorer()


default_sign_dict = {'Activation': 0,
                     'Inhibition': 1,
                     'IncreaseAmount': 0,
                     'DecreaseAmount': 1}

INDRA_ROOT = path.abspath(path.dirname(path.abspath(indra.__file__)))
INDRA_RESOURCES = path.join(INDRA_ROOT, 'resources')
with open(path.join(INDRA_RESOURCES, 'source_mapping.json'), 'r') as f:
    db_source_mapping = json.load(f)


class IndraNet(nx.MultiDiGraph):
    """A Networkx representation of INDRA Statements."""
    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self._is_multi = True
        self.mandatory_columns = ['agA_name', 'agB_name', 'agA_ns', 'agA_id',
                                  'agB_ns', 'agB_id', 'stmt_type',
                                  'evidence_count', 'stmt_hash', 'belief',
                                  'source_counts']

    @classmethod
    def from_df(cls, df):
        """Create an IndraNet MultiDiGraph from a pandas DataFrame.

        Returns an instance of IndraNet with graph data filled out from a
        dataframe containing pairwise interactions.

        Parameters
        ----------
        df : pd.DataFrame
            A :py:class:`pandas.DataFrame` with each row containing node and
            edge data for one edge. Indices are used to distinguish
            multiedges between a pair of nodes. Any columns not part of the
            below mentioned mandatory columns are considered extra attributes.
            Columns starting with 'agA\_' or 'agB\_' (excluding the agA/B_name)
            will be added to its respective nodes as node attributes. Any other
            columns will be added as edge attributes.

            Mandatory columns are : `agA_name`, `agB_name`, `agA_ns`, `agA_id`,
            `agB_ns`, `agB_id`, `stmt_type`, `evidence_count`, `stmt_hash`,
            `belief` and `source_counts`.

        Returns
        -------
        IndraNet
            An IndraNet object
        """
        graph = cls()
        mandatory_columns = graph.mandatory_columns
        if not set(mandatory_columns).issubset(set(df.columns)):
            raise ValueError('Missing one or more columns of %s in data '
                             'frame' % mandatory_columns)
        node_keys = {'agA': set(), 'agB': set()}
        edge_keys = set()
        for key in df.columns:
            if key not in mandatory_columns:
                if key.startswith('agA_'):
                    node_keys['agA'].add(key)
                if key.startswith('agB_'):
                    node_keys['agB'].add(key)
                if not key.startswith('ag'):
                    edge_keys.add(key)
        index = 0
        skipped = 0
        for index, row in df.iterrows():
            if row['agA_name'] is None or row['agB_name'] is None:
                skipped += 1
                logger.warning('None found as node (index %d)' % index)
                continue
            # Check and get node/edge attributes
            nodeA_attr = {}
            nodeB_attr = {}
            edge_attr = {}
            if node_keys['agA']:
                for key in node_keys['agA']:
                    nodeA_attr[key] = row[key]
            if node_keys['agB']:
                for key in node_keys['agB']:
                    nodeB_attr[key] = row[key]
            if edge_keys:
                for key in edge_keys:
                    edge_attr[key] = row[key]

            # Add non-existing nodes
            if row['agA_name'] not in graph.nodes:
                graph.add_node(row['agA_name'], ns=row['agA_ns'],
                               id=row['agA_id'], **nodeA_attr)
            if row['agB_name'] not in graph.nodes:
                graph.add_node(row['agB_name'], ns=row['agB_ns'],
                               id=row['agB_id'], **nodeB_attr)
            # Add edges
            ed = {'u_for_edge': row['agA_name'],
                  'v_for_edge': row['agB_name'],
                  'stmt_hash': row['stmt_hash'],
                  'stmt_type': row['stmt_type'],
                  'evidence_count': row['evidence_count'],
                  'belief': row['belief'],
                  'source_counts': row['source_counts'],
                  **edge_attr}
            graph.add_edge(**ed)
        if skipped:
            logger.warning('Skipped %d edges with None as node' % skipped)
        return graph

    def to_digraph(self, flattening_method=None, weight_mapping=None):
        """Flatten the IndraNet to a DiGraph

        Parameters
        ----------
        flattening_method : str|function
            The method to use when updating the belief for the flattened edge
        weight_mapping : function
            A function taking at least the graph G as an argument and
            returning G after adding edge weights as an edge attribute to the
            flattened edges using the reserved keyword 'weight'.

        Returns
        -------
        G : IndraNet(nx.DiGraph)
            An IndraNet graph flattened to a DiGraph
        """
        G = nx.DiGraph()
        for u, v, data in self.edges(data=True):
            # Add nodes and their attributes
            if u not in G.nodes:
                G.add_node(u, **self.nodes[u])
            if v not in G.nodes:
                G.add_node(v, **self.nodes[v])
            # Add edges and their attributes
            if G.has_edge(u, v):
                G[u][v]['statements'].append(data)
            else:
                G.add_edge(u, v, statements=[data])
        G = self._update_edge_belief(G, flattening_method)
        if weight_mapping:
            G = weight_mapping(G)
        return G

    def to_signed_graph(self, sign_dict=None,
                        flattening_method=None, weight_mapping=None):
        """Flatten the IndraNet to a signed graph.

        Parameters
        ----------
        sign_dict : dict
            A dictionary mapping a Statement type to a sign to be used for
            the edge. By default only Activation and IncreaseAmount are added
            as positive edges and Inhibition and DecreaseAmount are added as
            negative edges, but a user can pass any other Statement types in
            a dictionary.
        flattening_method : str or function(networkx.DiGraph, edge)
            The method to use when updating the belief for the flattened edge.

            If a string is provided, it must be one of the predefined options
            'simple_scorer' or 'complementary_belief'.

            If a function is provided, it must take the flattened graph 'G'
            and an edge 'edge' to perform the belief flattening on and return
            a number:

            >>> def flattening_function(G, edge):
            ...     # Return the average belief score of the constituent edges
            ...     all_beliefs = [s['belief']
            ...         for s in G.edges[edge]['statements']]
            ...     return sum(all_beliefs)/len(all_beliefs)

        weight_mapping : function(networkx.DiGraph)
            A function taking at least the graph G as an argument and
            returning G after adding edge weights as an edge attribute to the
            flattened edges using the reserved keyword 'weight'.

            Example:

            >>> def weight_mapping(G):
            ...     # Sets the flattened weight to the average of the
            ...     # inverse source count
            ...     for edge in G.edges:
            ...         w = [1/s['evidence_count']
            ...             for s in G.edges[edge]['statements']]
            ...         G.edges[edge]['weight'] = sum(w)/len(w)
            ...     return G

        Returns
        -------
        SG : IndraNet(nx.MultiDiGraph)
            An IndraNet graph flattened to a signed graph
        """
        sign_dict = default_sign_dict if not sign_dict else sign_dict

        SG = nx.MultiDiGraph()
        for u, v, data in self.edges(data=True):
            # Explicit 'is not None' needed to accept 0
            if data.get('initial_sign') is not None:
                sign = data['initial_sign']
            elif data['stmt_type'] not in sign_dict:
                continue
            else:
                sign = sign_dict[data['stmt_type']]
            if SG.has_edge(u, v, sign):
                SG[u][v][sign]['statements'].append(data)
            else:
                SG.add_edge(u, v, sign, statements=[data], sign=sign)
        SG = self._update_edge_belief(SG, flattening_method)
        if weight_mapping:
            SG = weight_mapping(SG)
        return SG

    @classmethod
    def digraph_from_df(cls, df, flattening_method=None, weight_mapping=None):
        """Create a digraph from a pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            The dataframe to build the graph from.
        flattening_method : str or function(networkx.DiGraph, edge)
            The method to use when updating the belief for the flattened edge.
        weight_mapping : function(networkx.DiGraph)
            A function taking at least the graph G as an argument and
            returning G after adding edge weights as an edge attribute to the
            flattened edges using the reserved keyword 'weight'.

        Returns
        -------
        IndraNet(nx.DiGraph)
             An IndraNet graph flattened to a DiGraph"""
        net = cls.from_df(df)
        return net.to_digraph(flattening_method=flattening_method,
                              weight_mapping=weight_mapping)

    @classmethod
    def signed_from_df(cls, df, sign_dict=None, flattening_method=None,
                       weight_mapping=None):
        """Create a signed graph from a pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            The dataframe to build the signed graph from.
        sign_dict : dict
            A dictionary mapping a Statement type to a sign to be used for
            the edge. By default only Activation and IncreaseAmount are added
            as positive edges and Inhibition and DecreaseAmount are added as
            negative edges, but a user can pass any other Statement types in
            a dictionary.
        flattening_method : str or function(networkx.DiGraph, edge)
            The method to use when updating the belief for the flattened edge.
        weight_mapping : function(networkx.DiGraph)
            A function taking at least the graph G as an argument and
            returning G after adding edge weights as an edge attribute to the
            flattened edges using the reserved keyword 'weight'.

        Returns
        -------
        IndraNet(nx.MultiDiGraph)
            An IndraNet graph flattened to a signed graph
        """
        net = cls.from_df(df)
        return net.to_signed_graph(sign_dict=sign_dict,
                                   flattening_method=flattening_method,
                                   weight_mapping=weight_mapping)

    @staticmethod
    def _update_edge_belief(G, flattening_method):
        """G must be or be a child of an nx.Graph object. If
        'flattening_method' is a function, it must take at least the graph G
        and an edge and return a number (the new belief for the flattened
        edge).

        We assume that G is the flattened graph and that all its edges have an
        edge attribute called 'statements' containing a list of dictionaries
        representing the edge data of all the edges in the un-flattened graph
        that were mapped to the corresponding flattened edge in G.
        """

        if not flattening_method or flattening_method == 'simple_scorer':
            for e in G.edges:
                G.edges[e]['belief'] = _simple_scorer_update(G, edge=e)
        elif flattening_method == 'complementary_belief':
            for e in G.edges:
                G.edges[e]['belief'] = _complementary_belief(G, edge=e)
        else:
            for e in G.edges:
                G.edges[e]['belief'] = flattening_method(G, edge=e)
        return G


def _simple_scorer_update(G, edge):
    evidence_list = []
    for stmt_data in G.edges[edge]['statements']:
        for k, v in stmt_data['source_counts'].items():
            if k in db_source_mapping:
                s = db_source_mapping[k]
            else:
                s = k
            for _ in range(v):
                evidence_list.append(Evidence(source_api=s))
    return simple_scorer.score_statement(st=Statement(evidence=evidence_list))


def _complementary_belief(G, edge):
    # Aggregate belief score: 1-prod(1-belief_i)
    np.seterr(all='raise')
    NP_PRECISION = 10 ** -np.finfo(np.longfloat).precision  # Numpy precision
    belief_list = [s['belief'] for s in G.edges[edge]['statements']]
    try:
        ag_belief = np.longfloat(1.0) - np.prod(np.fromiter(
            map(lambda belief: np.longfloat(1.0) - belief, belief_list),
            dtype=np.longfloat))
    except FloatingPointError as err:
        logger.warning('%s: Resetting ag_belief to 10*np.longfloat precision '
                       '(%.0e)' % (err, Decimal(NP_PRECISION * 10)))
        ag_belief = NP_PRECISION * 10
    return ag_belief
