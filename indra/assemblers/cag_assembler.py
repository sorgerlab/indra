# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function, unicode_literals
from builtins import object, dict, str
import os
import json
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
    def __init__(self, stmts=None, grounding_threshold=None):
        if not stmts:
            self.statements = []
        else:
            self.statements = stmts
        self.grounding_threshold = grounding_threshold

    def add_statements(self, stmts):
        """Add a list of Statements to the assembler."""
        self.statements += stmts

    def make_model(self, grounding_threshold=None):
        """Return a networkx MultiDiGraph representing a causal analysis graph.

        Parameters
        ----------
        grounding_threshold: Optional[float]
            Minimum threshold score for Eidos grounding.

        Returns
        -------
        nx.MultiDiGraph
            The assembled CAG.
        """
        self.grounding_threshold = grounding_threshold

        # Filter to Influence Statements which are currently supported
        statements = [stmt for stmt in self.statements if
                      isinstance(stmt, Influence)]

        # Initialize graph
        self.CAG = nx.MultiDiGraph()

        # Add nodes and edges to the graph
        for s in statements:
            # Get standardized name of subject and object
            # subj, obj = (self._node_name(s.subj), self._node_name(s.obj))

            # See if both subject and object have polarities given
            has_both_polarity = (s.subj_delta['polarity'] is not None and
                                 s.obj_delta['polarity'] is not None)

            # Add the nodes to the graph
            for node, delta in zip((s.subj, s.obj),
                                   (s.subj_delta, s.obj_delta)):
                self.CAG.add_node(self._node_name(node),
                        simulable=has_both_polarity, mods=delta['adjectives'])

            # Edge is solid if both nodes have polarity given
            linestyle = 'solid' if has_both_polarity else 'dotted'
            if has_both_polarity:
                same_polarity = s.subj_delta['polarity'] == s.obj_delta['polarity']
                if same_polarity:
                    targetArrowShape, linecolor = ('circle', 'green')
                else:
                    targetArrowShape, linecolor = ('tee', 'maroon')
            else:
                targetArrowShape, linecolor = ('triangle', 'maroon')


            # Add edge to the graph with metadata from statement
            provenance = []
            if s.evidence:
                provenance = s.evidence[0].annotations.get('provenance', [])
                if provenance:
                    provenance[0]['text'] = s.evidence[0].text
            self.CAG.add_edge(
                    self._node_name(s.subj),
                    self._node_name(s.obj),
                    subj_polarity    = s.subj_delta['polarity'],
                    subj_adjectives  = s.subj_delta['adjectives'],
                    obj_polarity     = s.obj_delta['polarity'],
                    obj_adjectives   = s.obj_delta['adjectives'],
                    linestyle        = linestyle,
                    linecolor=linecolor,
                    targetArrowShape=targetArrowShape,
                    provenance=provenance,
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
        def _create_edge_data_dict(e):
            """Return a dict from a MultiDiGraph edge for CytoscapeJS export."""
            # A hack to get rid of the redundant 'Provenance' label.
            if e[3].get('provenance'):
                tooltip = e[3]['provenance'][0]
                if tooltip.get('@type'):
                    del tooltip['@type']
            else:
                tooltip = None
            edge_data_dict = {
                    'id'               : e[0]+'_'+e[1],
                    'source'           : e[0],
                    'target'           : e[1],
                    'linestyle'        : e[3]["linestyle"],
                    'linecolor'        : e[3]["linecolor"],
                    'targetArrowShape' : e[3]["targetArrowShape"],
                    'subj_adjectives'  : e[3]["subj_adjectives"],
                    'subj_polarity'    : e[3]["subj_polarity"],
                    'obj_adjectives'   : e[3]["obj_adjectives"],
                    'obj_polarity'     : e[3]["obj_polarity"],
                    'tooltip'          : tooltip,
                    'simulable'        : False if (
                        e[3]['obj_polarity'] is None or
                        e[3]['subj_polarity'] is None) else True,
                    }
            return edge_data_dict

        return {
                'nodes': [{'data': {
                    'id': n[0],
                    'simulable': n[1]['simulable'],
                    'tooltip': 'Modifiers: '+json.dumps(n[1]['mods'])}
                    } for n in self.CAG.nodes(data=True)],

                'edges': [{'data': _create_edge_data_dict(e)}
                           for e in self.CAG.edges(data=True, keys=True)]
                }

    def generate_jupyter_js(self, cyjs_style=None, cyjs_layout=None):
        """Generate Javascript from a template to run in Jupyter notebooks.

        Parameters
        ----------
        cyjs_style : Optional[dict]
            A dict that sets CytoscapeJS style as specified in
            https://github.com/cytoscape/cytoscape.js/blob/master/documentation/md/style.md.

        cyjs_layout : Optional[dict]
            A dict that sets CytoscapeJS
            `layout parameters <http://js.cytoscape.org/#core/layout>`_.

        Returns
        -------
        str
            A Javascript string to be rendered in a Jupyter notebook cell.
        """
        # First, export the CAG to CyJS
        cyjs_elements = self.export_to_cytoscapejs()
        # Load the Javascript template
        tempf = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'cag_template.js')
        with open(tempf, 'r') as fh:
            template = fh.read()
        # Load the default style and layout
        stylef = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'cag_style.json')
        with open(stylef, 'r') as fh:
            style = json.load(fh)
        # Apply style and layout only if arg wasn't passed in
        if cyjs_style is None:
            cyjs_style = style['style']
        if cyjs_layout is None:
            cyjs_layout = style['layout']
        # Now fill in the template
        formatted_args = tuple(json.dumps(x, indent=2) for x in
                               (cyjs_elements, cyjs_style, cyjs_layout))
        js_str = template % formatted_args
        return js_str

    def _node_name(self, concept):
        """Return a standardized name for a node given a Concept."""
        if (# grounding threshold is specified
            self.grounding_threshold is not None
            # Eidos groundings are present
            and concept.db_refs['EIDOS']
            # The grounding score is above the grounding threshold
            and concept.db_refs['EIDOS'][0][1] > self.grounding_threshold):
                entry = concept.db_refs['EIDOS'][0][0]
                return entry.split('/')[-1].replace('_', ' ').capitalize()
        else:
            return concept.name.capitalize()
