from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, bytes
from past.builtins import basestring
import json
import logging
import requests
from .processor import EidosProcessor

logger = logging.getLogger(__name__)


try:
    # For text reading
    from .reader import EidosReader
    eidos_reader = EidosReader()
except Exception as e:
    logger.warning('Could not instantiate Eidos reader, text reading '
                   'will not be available.')
    eidos_reader = None


def process_text(text, save_json='eidos_output.json',
                 webservice=None, grounding_ns=None):
    """Return an EidosProcessor by processing the given text.

    This constructs a reader object via Java and extracts mentions
    from the text. It then serializes the mentions into JSON and
    processes the result with process_json.

    Parameters
    ----------
    text : str
        The text to be processed.
    save_json : Optional[str]
        The name of a file in which to dump the JSON output of Eidos.
    webservice : Optional[str]
        An Eidos reader web service URL to send the request to.
        If None, the reading is assumed to be done with the Eidos JAR rather
        than via a web service. Default: None
    grounding_ns : Optional[list]
        A list of name spaces for which INDRA should represent groundings, when
        given. If not specified or None, all grounding name spaces are
        propagated. If an empty list, no groundings are propagated.
        Example: ['UN', 'WM'], Default: None

    Returns
    -------
    ep : EidosProcessor
        An EidosProcessor containing the extracted INDRA Statements in its
        statements attribute.
    """
    if not webservice:
        if eidos_reader is None:
            logger.error('Eidos reader is not available.')
            return None
        json_dict = eidos_reader.process_text(text)
    else:
        if webservice.endswith('/'):
            webservice = webservice[:-1]
        res = requests.post('%s/process_text' % webservice,
                            json={'text': text})
        json_dict = res.json()
    if save_json:
        with open(save_json, 'wt') as fh:
            json.dump(json_dict, fh, indent=2)
    return process_json(json_dict, grounding_ns=grounding_ns)


def process_json_file(file_name, grounding_ns=None):
    """Return an EidosProcessor by processing the given Eidos JSON-LD file.

    This function is useful if the output from Eidos is saved as a file and
    needs to be processed.

    Parameters
    ----------
    file_name : str
        The name of the JSON-LD file to be processed.
    grounding_ns : Optional[list]
        A list of name spaces for which INDRA should represent groundings, when
        given. If not specified or None, all grounding name spaces are
        propagated. If an empty list, no groundings are propagated.
        Example: ['UN', 'WM'], Default: None

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in its statements attribute.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_str(json_str, grounding_ns=grounding_ns)
    except IOError:
        logger.exception('Could not read file %s.' % file_name)


def process_json_str(json_str, grounding_ns=None):
    """Return an EidosProcessor by processing the Eidos JSON-LD string.

    Parameters
    ----------
    json_str : str
        The JSON-LD string to be processed.
    grounding_ns : Optional[list]
        A list of name spaces for which INDRA should represent groundings, when
        given. If not specified or None, all grounding name spaces are
        propagated. If an empty list, no groundings are propagated.
        Example: ['UN', 'WM'], Default: None

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in its statements attribute.
    """
    json_dict = json.loads(json_str)
    return process_json(json_dict, grounding_ns=grounding_ns)


def process_json(json_dict, grounding_ns=None):
    """Return an EidosProcessor by processing a Eidos JSON-LD dict.

    Parameters
    ----------
    json_dict : dict
        The JSON-LD dict to be processed.
    grounding_ns : Optional[list]
        A list of name spaces for which INDRA should represent groundings, when
        given. If not specified or None, all grounding name spaces are
        propagated. If an empty list, no groundings are propagated.
        Example: ['UN', 'WM'], Default: None

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in its statements attribute.
    """
    ep = EidosProcessor(json_dict, grounding_ns=grounding_ns)
    ep.extract_causal_relations()
    ep.extract_correlations()
    ep.extract_events()
    return ep


def initialize_reader():
    """Instantiate an Eidos reader for fast subsequent reading."""
    eidos_reader.process_text('')
