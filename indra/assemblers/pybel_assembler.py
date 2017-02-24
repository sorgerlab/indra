from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.databases import hgnc_client
import logging
import pybel
import pybel.constants as pc
from pybel.parser.modifiers import PmodParser

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('pybel_assembler')

class PybelAssembler(object):

    def __init__(self, stmts=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts

        self.model = None


    def make_model(self, **kwargs):
        self.model = pybel.BELGraph(**kwargs)
        for stmt in self.statements:
            # Convert phosphorylation statements into
            # kin(p(agent)) => pmod(p(agent))
            if isinstance(stmt, Phosphorylation):
                # Define the enzyme node
                (enz_node, enz_attr) = _get_agent_node(stmt.enz)
                (sub_node, sub_attr) = _get_agent_node(stmt.sub)
                sub_attr[pc.VARIANTS] = [
                                    {pc.KIND: pc.PMOD,
                                     pc.IDENTIFIER: {
                                       pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                                       pc.NAME: 'Ph'},
                                     PmodParser.CODE: stmt.residue,
                                     PmodParser.POSITION: stmt.position }]
                self.model.add_node(enz_node, enz_attr)
                self.model.add_node(sub_node, sub_attr)
        #braf = ('Protein', 'HGNC', 'BRAF')
        #mek = ('Protein', 'HGNC', 'MAP2K1')
        #belgraph.add_node(braf, {'function': 'Protein', 'name':'BRAF',
        #                         'namespace':'HGNC'})
        #belgraph.add_node(mek, {'function': 'Protein', 'name':'MAP2K1',
        #                        'namespace':'HGNC'})
        edge = {pc.ANNOTATIONS: {},
                pc.CITATION: {
                    'authors': '',
                    'comments': '',
                    'date': '',
                    'name': 'Mol Cell Biol 2001 Apr 21(8) 2659-70',
                    'reference': '11283246',
                    'type': 'PubMed'},
                pc.EVIDENCE: 'Some evidence text.',
                pc.SUBJECT: {
                    pc.EFFECT: {pc.NAME: 'kin',
                                pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE},
                    pc.MODIFIER: pc.ACTIVITY},
                pc.RELATION: pc.INCREASES,
                pc.OBJECT: {}}
                #'object': {
                #    'effect': {'name': 'kin', 'namespace': 'bel'},
                #    'modifier': 'Activity'}}
        self.model.add_edge(enz_node, sub_node, attr_dict=edge)
        return self.model

def _get_agent_node(agent):
    if agent is None:
        return (None, None)
    # FIXME: For now, deal only with agents having HGNC grounding
    # FIXME: For now, don't look at bound/mod/mut/loc/activity conditions
    hgnc_id = agent.db_refs.get('HGNC')
    if hgnc_id is None:
        return (None, None)
    hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
    node = (pc.PROTEIN, 'HGNC', hgnc_name)
    node_attr = {pc.FUNCTION: pc.PROTEIN, pc.NAMESPACE: 'HGNC',
                 pc.NAME: hgnc_name}
    return (node, node_attr)

def _get_citation(evidence):
    pass

"""
structure of bel graph.
- node/tuple, maps to
  - dict with node/tuples as keys, each one mapped to
    - a dict, with integers (negative integers???) as keys, each one mapped to
      - a dict with edge information?

inner dict has fields such as
{'annotations': {},
'citation': {
    'authors': '',
    'comments': '',
    'date': '',
    'name': 'Mol Cell Biol 2001 Apr 21(8) 2659-70',
    'reference': '11283246',
    'type': 'PubMed'},
'evidence': 'Modified assertion',
'subject': {
    'effect': {'name': 'kin', 'namespace': 'bel'},
    'modifier': 'Activity'}}}
'relation': 'increases',
'object': {
    'effect': {'name': 'kin', 'namespace': 'bel'},
    'modifier': 'Activity'},

Note that Activity modifiers are listed in the subject/object dicts.

Edges are added between top level nodes and modified nodes using a hasVariant
relation:

('Protein',
  'HGNC',
  'BRAF',
  ('hgvs', 'p.Val600Glu')): {-4: {'relation': 'hasVariant'}},


"""

# Need function for building up terms from agent conditions, including
# (mainly for now) modification conditions

# Example nodes in graph:
"""
('Protein', 'RGD', 'Pdpk1', ('pmod', ('bel', 'Ph'))),
 ('Complex', ('Protein', 'HGNC', 'JUP'), ('Protein', 'PFH', 'AXIN Family')),
 ('RNA', 'HGNC', 'LTV1'),
 ('RNA', 'MGI', 'Nfix'),
 ('Protein', 'HGNC', 'PAX7'),
 ('RNA', 'HGNC', 'NEUROG1'),
 ('RNA', 'MGI', 'Mup1'),
 ('Protein', 'HGNC', 'PPARA'),
 ('RNA', 'EGID', '66935'),
 ('RNA', 'HGNC', 'AKIRIN1'),
"""

# Then need patterns for assembling modifications, activations, rasgef/gap

# Can start with this, then add evidence/context, etc.


