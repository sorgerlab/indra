# -*- coding: utf-8 -*-

"""A processor for the `Drug Gene Interaction Database (DGI-DB) <http://www.dgidb.org>`_.

* `Integration of the Drug–Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts
   <https://doi.org/10.1093/nar/gkaa1084>`_. Freshour, *et al*. Nucleic Acids Research. 2020 Nov 25.

Interactions data from the January 2021 release can be obtained at the
following URLs:

* https://www.dgidb.org/data/monthly_tsvs/2021-Jan/interactions.tsv
"""

from .api import get_version_df, process_df, process_version
from .processor import DGIProcessor
