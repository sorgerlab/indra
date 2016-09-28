from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra import bel
from indra.util import unicode_strs

def test_bel_ndex_query():
    bp = bel.process_ndex_neighborhood(['NFKB1'])
    unicode_strs(bp.statements)

