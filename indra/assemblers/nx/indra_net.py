import networkx as nx
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class IndraNet(nx.MultiDiGraph):
    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self._is_multi = True

    def _add_node(self, name, ns, id, **attr):
        all_attr = {'ns': ns, 'id': id, **attr}
        self.add_node(node_for_adding=name, **all_attr)

    def _add_edge(self, u, v, key=None, **attr):
        self.add_edge(u_for_edge=u, v_for_edge=v, key=key, **attr)

    @classmethod
    def from_df(cls, df=pd.DataFrame(), belief_dict=None, strat_ev_dict=None,
                multi=True):
        """idea: enforce pair of ('agA_<attribute>', 'agB_<attribute>') for
        node attributes, otherwise considered edge attribute

        :param df: pd.DataFrame
        :param belief_dict:
        :param strat_ev_dict:
        :param multi:
        :return:
        """
        if not {'agA_name', 'agB_name', 'agA_ns', 'agA_id', 'agB_ns', 'agB_id',
                'stmt_type', 'evidence_count', 'hash', 'belief',
                'evidence'}.issubset(
             set(df.columns)):
            raise ValueError('Missing required columns in data frame')

        index = 0
        skipped = 0
        for index, row in df.iterrows():
            if row['agA_name'] is None or row['agB_name'] is None:
                skipped += 1
                logger.warning('None found as node (index %d)' % index)
                continue
            pass
        return cls

    @classmethod
    def to_type_graph(cls):
        # Will wrap 'from_df' and collapse edges as necessary
        cls.from_df()
        return cls

    def is_multigraph(self):
        return self._is_multi
