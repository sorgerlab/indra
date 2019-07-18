import networkx as nx
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class IndraNet(nx.MultiDiGraph):
    """IndraNet network object

    Methods
    -------
    from_df(df): (classmethod)
        Return an instance of IndraNet, with graph data filled out from a
        dataframe (df) containing pairwise interactions.
    """
    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self._is_multi = True

    @classmethod
    def from_df(cls, df):
        graph = cls()
        mandatory_columns = ['agA_name', 'agB_name', 'agA_ns', 'agA_id',
                             'agB_ns', 'agB_id', 'stmt_type', 'evidence_count',
                             'hash', 'belief']
        """Create an IndraNet MultiDiGraph from a pandas DataFrame

        Parameters
        ----------
        df: pd.DataFrame
            A data frame with each row containing node and edge data for one
            edge. Mandatory columns are {m}. Hashes are used to distinguish
            multiedges between a pair of nodes. Any other columns are
            considered extra node or edge attributes. Any columns starting
            with 'agA_' or 'agB_' (excluding the mandatory agA/B_name) will
            be added to its respective nodes as node attributes. Columns not
            starting with 'agA_' or 'agB_' will be added as edge attributes.

        Returns
        -------
        in : IndraNet
            An IndraNet object
        """.format(m=mandatory_columns)
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
                  'key': row['hash'],
                  'stmt_type': row['stmt_type'],
                  'evidence_count': row['evidence_count'],
                  'belief': row['belief'],
                  **edge_attr}
            graph.add_edge(**ed)
        if skipped:
            logger.warning('Skipped %d edges with None as node' % skipped)
        return graph
