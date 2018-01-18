from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, bytes
from past.builtins import basestring
import json
import logging
from .processor import EidosProcessor

logger = logging.getLogger('eidos')


try:
    # For text reading
    from .eidos_reader import EidosReader
    eidos_reader = EidosReader()
except Exception:
    logger.error('Could not instantiate Eidos reader, text reading '
                 'will not be available.')
    eidos_reader = None


def process_text(text, save_json='eidos_output.json'):
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

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    if eidos_reader is None:
        logger.error('Eidos reader is not available.')
        return None
    json_dict = eidos_reader.process_text(text)
    if save_json:
        with open(save_json, 'wt') as fh:
            json.dump(json_dict, fh, indent=1)
    return process_json(json_dict)


def process_json_file(file_name):
    """Return an EidosProcessor by processing the given Eidos json file.

    The output from the Eidos reader is in json format. This function is
    useful if the output is saved as a file and needs to be processed.

    Parameters
    ----------
    file_name : str
        The name of the json file to be processed.

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_str(json_str)
    except IOError:
        logger.error('Could not read file %s.' % file_name)


def process_json_str(json_str):
    """Return an EidosProcessor by processing the given Eidos json string.

    The output from the Eidos parser is in json format.

    Parameters
    ----------
    json_str : str
        The json string to be processed.

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    if not isinstance(json_str, basestring):
        raise TypeError('{} is {} instead of {}'.format(json_str,
                                                        json_str.__class__,
                                                        basestring))
    try:
        json_dict = json.loads(json_str)
    except ValueError:
        logger.error('Could not decode JSON string.')
        return None
    return process_json(json_dict)


def process_json(json_dict):
    """Return an EidosProcessor by processing the given Eidos json dict.

    Parameters
    ----------
    json_dict : dict
        The json dict to be processed.

    Returns
    -------
    ep : EidosProcessor
        A EidosProcessor containing the extracted INDRA Statements
        in ep.statements.
    """

    ep = EidosProcessor(json_dict)
    ep.get_events()
    return ep

