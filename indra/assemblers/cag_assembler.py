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
    def __init__(self, stmts=None):
        if not stmts:
            self.statements = []
        else:
            self.statements = stmts
        self.grounding_threshold = None
        self.grounding_ontology = 'UN'
        self.CAG = None

    def add_statements(self, stmts):
        """Add a list of Statements to the assembler."""
        self.statements += stmts

    def make_model(self, grounding_ontology='UN', grounding_threshold=None):
        """Return a networkx MultiDiGraph representing a causal analysis graph.

        Parameters
        ----------
        grounding_ontology : Optional[str]
            The ontology from which the grounding should be taken
            (e.g. UN, FAO)
        grounding_threshold : Optional[float]
            Minimum threshold score for Eidos grounding.

        Returns
        -------
        nx.MultiDiGraph
            The assembled CAG.
        """
        if grounding_threshold is not None:
            self.grounding_threshold = grounding_threshold

        self.grounding_ontology = grounding_ontology

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
                                  simulable=has_both_polarity,
                                  mods=delta['adjectives'])

            # Edge is solid if both nodes have polarity given
            linestyle = 'solid' if has_both_polarity else 'dotted'
            if has_both_polarity:
                same_polarity = (s.subj_delta['polarity'] ==
                                 s.obj_delta['polarity'])
                if same_polarity:
                    target_arrow_shape, linecolor = ('circle', 'green')
                else:
                    target_arrow_shape, linecolor = ('tee', 'maroon')
            else:
                target_arrow_shape, linecolor = ('triangle', 'maroon')

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
                    targetArrowShape=target_arrow_shape,
                    provenance=provenance,
                )

        return self.CAG

    def print_tsv(self, file_name):
        def _get_factor(stmt, concept, delta, evidence, raw_name):
            if evidence.source_api == 'eidos':
                if concept.db_refs[self.grounding_ontology]:
                    factor_norm = concept.db_refs[self.grounding_ontology][0][0]
                else:
                    factor_norm = ''
            elif evidence.source_api == 'hume':
                factor_norm = concept.db_refs['HUME'][0][0]
            elif evidence.source_api == 'cwms':
                factor_norm = concept.db_refs['CWMS']
            elif evidence.source_api == 'sofia':
                # TODO extract ontology catgory here
                factor_norm = concept.name
            mods = ', '.join(delta.get('adjectives', []))
            if delta.get('polarity') == -1:
                pol = 'decrease'
            elif delta.get('polarity') == 1:
                pol = 'increase'
            else:
                pol = ''
            name = raw_name if raw_name else concept.name
            return name, factor_norm, mods, pol

        def _get_evidence(evidence):
            # TODO: add sentence ID
            sent_id = ''
            location = evidence.annotations.get('Location')
            location = location if location is not None else ''
            time = evidence.annotations.get('Time')
            time = time if time is not None else ''
            ref = evidence.pmid if evidence.pmid is not None else ''
            return ref, evidence.source_api, sent_id, location, \
                time, evidence.text

        header = ['Source', 'System', 'Sentence ID',
                  'Factor A Text', 'Factor A Normalization',
                  'Factor A Modifiers', 'Factor A Polarity',
                  'Relation Text', 'Relation Normalization',
                  'Relation Modifiers',
                  'Factor B Text', 'Factor B Normalization',
                  'Factor B Modifiers', 'Factor B Polarity',
                  'Location', 'Time', 'Evidence',
                  'Relation ID']

        fh = open(file_name, 'w')
        fh.write('\t'.join(header) + '\n')

        # Filter to Influence Statements which are currently supported
        statements = [stmt for stmt in self.statements if
                      isinstance(stmt, Influence)]
        all_rows = []
        for idx, stmt in enumerate(statements):
            for evidence in stmt.evidence:
                source, system, sent_id, location, time, text = \
                    _get_evidence(evidence)
                factor_a, factor_a_norm, mod_a, pol_a = \
                    _get_factor(stmt, stmt.subj, stmt.subj_delta, evidence,
                                evidence.annotations['subj_text'])
                factor_b, factor_b_norm, mod_b, pol_b = \
                    _get_factor(stmt, stmt.obj, stmt.obj_delta, evidence,
                                evidence.annotations['obj_text'])
                relation_text = 'influences'
                # Can we get a more specific relation type here?
                relation_norm = ''
                relation_mod = ''
                row = [source, system, sent_id,
                       factor_a, factor_a_norm, mod_a, pol_a,
                       relation_text, relation_norm, relation_mod,
                       factor_b, factor_b_norm, mod_b, pol_b,
                       location, time, text, str(idx)]
                if row not in all_rows:
                    all_rows.append(row)
        for row in sorted(all_rows, key=lambda x: x[0]):
            fh.write('\t'.join(row) + '\n')
        fh.close()

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
            # The particular eidos ontology grounding (un/wdi/fao) is present
            and concept.db_refs[self.grounding_ontology]
            # The grounding score is above the grounding threshold
            and (concept.db_refs[self.grounding_ontology][0][1] >
                 self.grounding_threshold)):
                entry = concept.db_refs[self.grounding_ontology][0][0]
                return entry.split('/')[-1].replace('_', ' ').capitalize()
        else:
            return concept.name.capitalize()
