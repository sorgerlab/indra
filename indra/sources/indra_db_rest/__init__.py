"""The INDRA database client allows querying a web service that serves
content from a database of INDRA Statements collected and pre-assembled
from various sources.

Access to the webservice requires a URL (`INDRA_DB_REST_URL`) and
possibly an API key (`INDRA_DB_REST_API_KEY`), both of which may be placed in
your config file or as environment variables. If you do not have these
but would like to access the database REST API, you may contact the developers
to request a URL and API key.
"""
from .client_api import *
