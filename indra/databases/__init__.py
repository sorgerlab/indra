"""This module implements a number of clients for accessing and using
resources from biomedical entity databases and other third-party web
services that INDRA uses. Many of the resources these clients use are loaded
from resource files in the indra.resources module, in many cases also providing
access to web service endpoints."""

from .identifiers import get_identifiers_url, parse_identifiers_url, \
    url_prefixes
