__all__ = ['process_json_file_old', 'process_jsonld_file',
           'process_jsonld']

import json
import rdflib
import logging
from indra.sources.hume import processor
from os.path import abspath

logger = logging.getLogger('hume')


def process_jsonld_file(fname):
    """Process a JSON-LD file in the new format to extract Statements.

    Parameters
    ----------
    fname : str
        The path to the JSON-LD file to be processed.

    Returns
    -------
    indra.sources.hume.HumeProcessor
        A HumeProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    with open(fname, 'r') as fh:
        json_dict = json.load(fh)
    return process_jsonld(json_dict)


def process_jsonld(jsonld):
    """Process a JSON-LD string in the new format to extract Statements.

    Parameters
    ----------
    jsonld : dict
        The JSON-LD object to be processed.

    Returns
    -------
    indra.sources.hume.HumeProcessor
        A HumeProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    hp = processor.HumeJsonLdProcessor(jsonld)
    hp.get_events()
    return hp


def process_json_file_old(fname):
    """Process a JSON-LD file in the old format to extract Statements.

    Parameters
    ----------
    fname : str
        The path to the JSON-LD file to be processed.

    Returns
    -------
    bp : indra.sources.hume.HumeProcessor
        A HumeProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    graph = _load_graph(fname)
    bp = processor.HumeProcessor(graph)
    bp.get_statements()
    return bp


def _load_graph(fname):
    fname = abspath(fname)
    g = rdflib.Graph()
    with open(fname, 'rb') as fh:
        logger.info('Started loading graph from %s' % fname)
        g.parse(fh, format='json-ld')
        logger.info('Finished loading graph')
    return g
