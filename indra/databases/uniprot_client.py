from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging
from indra.util import read_unicode_csv

logger = logging.getLogger(__name__)


from protmapper.uniprot_client import *


def _build_uniprot_subcell_loc():
    fname = os.path.dirname(os.path.abspath(__file__)) +\
                '/../resources/uniprot_subcell_loc.tsv'
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    # Skip the header row
    up_to_go = {}
    for row in csv_rows:
        upid = row[0]
        goid = row[1]
        up_to_go[upid] = goid
    return up_to_go


uniprot_subcell_loc = _build_uniprot_subcell_loc()
