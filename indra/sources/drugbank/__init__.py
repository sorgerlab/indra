"""This module provides an API and processor for DrugBank content.
It builds on the XML-formatted data schema of DrugBank and expects
the XML file to be available locally. The full DrugBank download
can be obtained at: https://www.drugbank.ca/releases/latest. Once the
XML file is decompressed, it can be processed using the process_xml
function.

Alternatively, the latest DrugBank data can be automatically loaded via
:mod:`drugbank_downloader` with the following code after doing
``pip install drugbank_downloader bioversions``:

.. code-block:: python

    import pickle
    import indra.sources.drugbank
    processor = indra.sources.drugbank.get_drugbank_processor()
    with open('drugbank_indra_statements.pkl', 'wb') as file:
        pickle.dump(processor.statements, file, protocol=pickle.HIGHEST_PROTOCOL)
"""
from .api import get_drugbank_processor, process_xml
