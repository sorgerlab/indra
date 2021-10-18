# -*- coding: utf-8 -*-

"""Processor for the `CRowd Extracted Expression of Differential Signatures (CREEDS)
<https://maayanlab.cloud/CREEDS>`_ from [wang2016]_.

Contains results of differential gene expression experiments extracted
from the Gene Expression Omnibus due to three types of perturbations:

1. Single gene knockdown/knockout
2. Drug
3. Disease

.. [wang2016] Wang, Z., *et al*. (2016). `Extraction and analysis of signatures from the Gene
   Expression Omnibus by the crowd <https://doi.org/10.1038/ncomms12846>`_.
   *Nature Communications*, **7**(1), 12846.
"""

from .api import process_from_file, process_from_web
from .processor import (
    CREEDSChemicalProcessor,
    CREEDSDiseaseProcessor,
    CREEDSGeneProcessor,
)
