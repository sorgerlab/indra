# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, unicode_literals
from builtins import object, dict, str
import logging
import networkx as nx
from indra.statements import Influence

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cag_assembler')


class CAGAssembler(object):
    """Assembles a causal analysis graph from INDRA Statements.

    Parameters
    ----------
    stmts : Optional[list[indra.statement.Statements]]
        A list of INDRA Statements to be assembled. Currently supports
        Influence Statements.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    CAG : nx.MultiDiGraph
        A networkx MultiDiGraph object representing the causal analysis graph.
    """
    def __init__(self, stmts=None):
        if not stmts:
            self.statements = []
        else:
            self.statements = stmts

    def add_statements(self, stmts):
        """Add a list of Statements to the assembler."""
        self.statements += stmts

    def make_model(self):
        """Return a networkx MultiDiGraph representing a causal analysis graph.

        Returns
        -------
        nx.MultiDiGraph
            The assembled CAG.
        """
        # Filter to Influence Statements which are currently supported
        statements = [stmt for stmt in self.statements if
                      isinstance(stmt, Influence)]

        # Initialize graph
        self.CAG = nx.MultiDiGraph()

        # Extract node names from the INDRA Statements and add to graph
        node_names = set()
        for stmt in self.statements:
            for agent in stmt.agent_list():
                node_name = self._node_name(agent.name)
                if node_name not in self.CAG:
                    self.CAG.add_node(node_name, simulable=False)

        # Add edges to the graph
        for s in statements:
            subj, obj = (self._node_name(s.subj.name),
                         self._node_name(s.obj.name))

            if (s.subj_delta['polarity'] is not None and
                    s.obj_delta['polarity'] is not None):
                self.CAG.node[subj]['simulable'] = True
                self.CAG.node[obj]['simulable'] = True

            linestyle = 'solid' if (s.subj_delta['polarity'] is not None and
                                    s.obj_delta['polarity'] is not None) \
                        else 'dotted'

            self.CAG.add_edge(subj, obj, 
                    subj_polarity   = s.subj_delta['polarity'],
                    subj_adjectives = s.subj_delta['adjectives'],
                    obj_polarity    = s.obj_delta['polarity'],
                    obj_adjectives  = s.obj_delta['adjectives'],
                    linestyle       = linestyle
                    )

        return self.CAG

    def export_to_cytoscapejs(self):
        """Return CAG in format readable by CytoscapeJS.

        Return
        ------
        dict
            A JSON-like dict representing the graph for use with
            CytoscapeJS.
        """
        return {
                'nodes': [{'data': {'id': n[0],
                                    'simulable': n[1]['simulable']}}
                          for n in self.CAG.nodes(data=True)],
                'edges': [{
                    'data': {
                        'id'              : e[0]+'_'+e[1],
                        'source'          : e[0],
                        'target'          : e[1],
                        'linestyle'       : e[3]["linestyle"],
                        'subj_adjectives' : e[3]["subj_adjectives"],
                        'subj_polarity'   : e[3]["subj_polarity"],
                        'obj_adjectives'  : e[3]["obj_adjectives"],
                        'obj_polarity'    : e[3]["obj_polarity"],
                        'simulable' : False if (
                            e[3]['obj_polarity'] is None or
                            e[3]['subj_polarity'] is None) else True
                        }
                    }
                    for e in self.CAG.edges(data=True, keys=True)
                    ]
                }

    @staticmethod
    def _node_name(agent_name):
        """Return a standardized name for a node given an Agent name."""
        return agent_name.capitalize()
