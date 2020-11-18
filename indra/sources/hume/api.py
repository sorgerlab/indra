__all__ = ['process_jsonld_file', 'process_jsonld']

import json
import logging
from indra.sources.hume import processor

logger = logging.getLogger(__name__)

default_grounding_mode = 'flat'


def process_jsonld_file(fname, extract_filter=None,
                        grounding_mode=default_grounding_mode):
    """Process a JSON-LD file in the new format to extract Statements.

    Parameters
    ----------
    fname : str
        The path to the JSON-LD file to be processed.
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None
    grounding_mode : Optional[str]
        Selects whether 'flat' or 'compositional' groundings should be
        extracted. Default: 'flat'.

    Returns
    -------
    indra.sources.hume.HumeProcessor
        A HumeProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    with open(fname, 'r', encoding='utf-8') as fh:
        json_dict = json.load(fh)
    return process_jsonld(json_dict, extract_filter=extract_filter,
                          grounding_mode=grounding_mode)


def process_jsonld(jsonld, extract_filter=None,
                   grounding_mode=default_grounding_mode):
    """Process a JSON-LD string in the new format to extract Statements.

    Parameters
    ----------
    jsonld : dict
        The JSON-LD object to be processed.
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None
    grounding_mode : Optional[str]
        Selects whether 'flat' or 'compositional' groundings should be
        extracted. Default: 'flat'.

    Returns
    -------
    indra.sources.hume.HumeProcessor
        A HumeProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    if grounding_mode == 'flat':
        hp = processor.HumeJsonLdProcessor(jsonld)
    elif grounding_mode == 'compositional':
        hp = processor.HumeJsonLdProcessorCompositional(jsonld)
    else:
        raise ValueError('Invalid grounding mode: %s' % grounding_mode)

    if extract_filter is None or 'influence' in extract_filter:
        hp.extract_relations()
    if extract_filter is None or 'event' in extract_filter:
        hp.extract_events()
    return hp
