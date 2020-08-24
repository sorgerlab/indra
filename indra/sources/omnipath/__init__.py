"""
The OmniPath module accesses biomolecular interaction data from various
curated databases using the OmniPath API (see
https://saezlab.github.io/pypath/html/index.html#webservice) and processes
the returned data into statements using the OmniPathProcessor.

Currently, the following data is collected:
  - Modifications from the PTMS endpoint https://saezlab.github.io/pypath/html/index.html#enzyme-substrate-interactions
  - Ligand-Receptor data from the interactions endpoint https://saezlab.github.io/pypath/html/index.html#interaction-datasets

To process all statements, use the function `process_from_web`:

>>> from indra.sources.omnipath import process_from_web
>>> omnipath_processor = process_from_web()
>>> stmts = omnipath_processor.statements
"""
from .api import process_from_web
from .processor import OmniPathProcessor
