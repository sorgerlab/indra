from __future__ import unicode_literals
from builtins import dict, str

from indra.databases import omnipath as op

def test_query_ptms():
    op.get_ptms(['MAP2K1'])

