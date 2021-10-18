# -*- coding: utf-8 -*-

"""Processor for the `CREEDS
<https://maayanlab.cloud/CREEDS>`_.

Contains results of differential gene expression experiments extracted
from the Gene Expression Omnibus due to three types of perturbations:

1. Single gene knockdown/knockout
2. Drug
3. Disease
"""

from .api import process_from_web
from .processor import CREEDSChemicalProcessor, CREEDSDiseaseProcessor, CREEDSGeneProcessor
