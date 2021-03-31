"""
This module implements getting content from the Human Protein Reference
Database (HPRD), a curated protein data resource, as INDRA Statements.
In particular, the module supports extracting post-translational
modifications, protein complexes, and (binary) protein-protein interactions
from HPRD.

More information about HPRD can be obtained at
`http://www.hprd.org <http://www.hprd.org>`_ and in these publications:

* Peri, S. et al. (2003). Development of Human Protein Reference Database as an
  initial platform for approaching systems biology in humans. Genome Research.
  13, 2363-2371.
* Prasad, T. S. K. et al. (2009). Human Protein Reference Database - 2009
  Update. Nucleic Acids Research. 37, D767-72.

Data from the final release of HPRD (version 9) can be obtained at the
following URLs:

* http://www.hprd.org/RELEASE9/HPRD_FLAT_FILES_041310.tar.gz (text files)
* http://www.hprd.org/RELEASE9/HPRD_XML_041310.tar.gz (XML)

This module is designed to process the text files obtained from the first
link listed above.
"""
from .processor import HprdProcessor
from .api import process_flat_files
