"""This module implements an interface to Biofactoid (https://biofactoid.org/)
which contains interactions curated from publications by
authors. Documents are retrieved from the web and processed into
INDRA Statements."""

from .api import process_json, process_from_web

