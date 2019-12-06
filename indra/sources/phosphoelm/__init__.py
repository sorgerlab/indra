"""This module provides an interface to the Phospho.ELM database and extracts
phosphorylation relationships as INDRA Statements. Phospho.ELM is available at
http://phospho.elm.eu.org/, see also
https://academic.oup.com/nar/article/39/suppl_1/D261/2506728"""
from .api import process_from_dump
