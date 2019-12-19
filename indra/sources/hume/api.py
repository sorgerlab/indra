__all__ = ['process_jsonld_file', 'process_jsonld']

import json
import logging
from indra.sources.hume import processor

logger = logging.getLogger(__name__)


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
    with open(fname, 'r', encoding='utf-8') as fh:
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
    hp.extract_relations()
    hp.extract_events()
    return hp
