from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import pybel

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

    def make_model(self):
        belgraph = pybel.BELGraph()
        braf = ('Protein', 'HGNC', 'BRAF')
        mek = ('Protein', 'HGNC', 'MAP2K1')
        edge = {'annotations': {},
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
                    'modifier': 'Activity'},
                'relation': 'increases',
                'object': {
                    'effect': {'name': 'kin', 'namespace': 'bel'},
                    'modifier': 'Activity'}}
        belgraph.add_edges_from([(braf, mek, edge)])
        self.model = belgraph
        return self.model

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


