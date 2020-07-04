"""This module provides an API and processor for DrugBank content.
It builds on the XML-formatted data schema of DrugBank and expects
the XML file to be available locally. The full DrugBank download
can be obtained at: https://www.drugbank.ca/releases/latest. Once the
XML file is decompressed, it can be processed using the process_xml
function."""
from .api import process_xml