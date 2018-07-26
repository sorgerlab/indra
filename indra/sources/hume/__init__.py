"""Hume is a general purpose reading system developed by BBN.

Currently, INDRA can process JSON-LD files produced by Hume. When available,
the API will be extended with access to the reader as a service.
"""

from .api import process_json_file_old, process_jsonld_file
