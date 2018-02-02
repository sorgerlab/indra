from __future__ import absolute_import, print_function, unicode_literals
from past.builtins import basestring
from builtins import object, dict, str
import json
import logging
import itertools
import numpy as np
from indra.statements import *
from indra.preassembler import Preassembler
import networkx as nx
from indra.util import flatMap
from indra.statements import Influence

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cyjs_assembler')


class CAGAssembler(object):
    """This class assembles a causal analysis graph from a set of INDRA 
    statements.

    Parameters
    ----------
    statements : Optional[List[Influence]]
        A list of INDRA Influence Statements to be assembled.

    Attributes
    ----------
    statements : List[Influence]
        A list of INDRA Statements to be assembled.
    CAG : nx.MultiDiGraph
        A networkx MultiDiGraph object representing the causal analysis graph.
    """
    def __init__(self, stmts):
        self.statements = stmts
        self.CAG = self.create_causal_analysis_graph()

    def create_causal_analysis_graph(self):
        """ Create a networkx MultiDiGraph object that represents a causal
        analysis graph. 
        
        Returns
        -------
        nx.MultiDiGraph

        """

        # Extract unique factors from the INDRA statements
        factors = set(flatMap(lambda x: (x.subj.name, x.obj.name),
                              self.statements))

        # Interleave partial deriviatives w.r.t. time to create index of the
        # latent state components.
        s_index = flatMap(lambda a: (a, '∂('+a+')/∂t'), sorted(factors))

        G = nx.MultiDiGraph() 

        for latent_state_component in s_index:
            G.add_node(latent_state_component.capitalize(), simulable=False)

        for s in self.statements:
            subj, obj = s.subj.name.capitalize(), s.obj.name.capitalize()

            if s.subj_delta['polarity'] != None and s.obj_delta['polarity'] != None:
                G.nodes[subj]['simulable'] = True
                G.nodes[obj]['simulable'] = True

            key = G.add_edge(subj, obj,
                        subj_polarity = s.subj_delta['polarity'],
                        subj_adjectives = s.subj_delta['adjectives'],
                        obj_polarity = s.obj_delta['polarity'],
                        obj_adjectives = s.obj_delta['adjectives'],
                        linestyle='dotted'
                    )
            if s.subj_delta['polarity'] != None and s.obj_delta['polarity'] != None:
                G[subj][obj][key]['linestyle']='solid'

        return G

    def export_to_cytoscapejs(self):
        """ Export networkx to format readable by CytoscapeJS """
        return {
                'nodes':[{'data':{'id': n[0], 'simulable': n[1]['simulable']}} 
                    for n in self.G.nodes(data=True)],
                'edges':[ {
                    'data': {
                        'id'              : e[0]+'_'+e[1],
                        'source'          : e[0],
                        'target'          : e[1],
                        'linestyle'       : e[3]["linestyle"],
                        'subj_adjectives' : e[3]["subj_adjectives"],
                        'subj_polarity'   : e[3]["subj_polarity"],
                        'obj_adjectives'  : e[3]["obj_adjectives"],
                        'obj_polarity'    : e[3]["obj_polarity"],
                        'simulable' : False if (e[3]['obj_polarity'] == None or
                            e[3]['subj_polarity'] == None) else True
                        }
                    } 
                    for e in self.G.edges(data=True, keys = True)
                    ]
                }


