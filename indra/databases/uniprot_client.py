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
    try:
        csv_rows = read_unicode_csv(fname, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        subcell_loc = {}
        for row in csv_rows:
            loc_id = row[0]
            loc_alias = row[3]
            subcell_loc[loc_id] = loc_alias
    except IOError:
        subcell_loc = {}
    return subcell_loc


uniprot_subcell_loc = _build_uniprot_subcell_loc()
