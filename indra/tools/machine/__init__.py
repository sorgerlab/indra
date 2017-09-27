"""
Prerequisites
-------------
First, install the machine-specific dependencies:

.. code-block:: sh

    pip install indra[machine]

Starting a New Model
--------------------
To start a new model, run

.. code-block:: sh

    python -m indra.tools.machine make model_name

Alternatively, the command line interface can be invoked with

.. code-block:: sh

    indra-machine make model_name

where model_name corresponds to the name of the model to initialize.

This script generates the following folders and files

- ``model_name``
- ``model_name/log.txt``
- ``model_name/config.yaml``
- ``model_name/jsons/``

You should the edit model\_name/config.yaml to set up the search terms
and optionally the credentials to use Twitter, Gmail or NDEx bindings.

Setting Up Search Terms
-----------------------
The config.yml file is a standard YAML configuration file. A template
is available in  model\_name/config.yaml after having created the machine.

Two important fields in config.yml are ``search_terms`` and ``search_genes``
both of which are YAML lists. The entries of ``search_terms`` are used
_directly_ as queries in PubMed search (for more information on PubMed
search strings,
read https://www.ncbi.nlm.nih.gov/books/NBK3827/#pubmedhelp.Searching\_PubMed).

Example:

.. code-block:: yaml

    search_terms:
    - breast cancer
    - proteasome
    - apoptosis

The entries of ``search_genes`` is a special list in which _only_ standard
HGNC gene symbols are allowed. Entries in this list are also used
to search PubMed but also serve as a list of `prior` genes that are known
to be relevant for the model.

#Entries in this can be used to search
#PubMed specifically for articles that are tagged with the gene's unique
#identifier rather than its string name. This mode of searching for articles
#on specific genes is much more reliable than searching for them using
#string names.

Example:

.. code-block:: yaml

    search_genes:
    - AKT1
    - MAPK3
    - EGFR

Extending a Model
-----------------
To extend a model, run

.. code-block:: sh

    python -m indra.tools.machine run_with_search model_name

Alternatively, the command line interface can be invoked with

.. code-block:: sh

    indra-machine run_with_search model_name

Extending a model involves extracting PMIDs from emails (if Gmail credentials
are given), and searching using INDRA's PubMed client with each entry
of search\_terms in config.yaml as a search term.  INDRA's literature
client is then used to find the full text corresponding to each PMID
or its abstract when the full text is not available.
The REACH parser is then used to read each new paper.
INDRA uses the REACH output to construct Statements corresponding to
mechanisms. It then adds them to an incremental model through a process of
assembly involving duplication and overlap resolution and the application of
filters.
"""

from . import config
from .config import *

__all__ = (
    config.__all__
)
