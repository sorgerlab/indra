from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, bytes
from past.builtins import basestring
import json
import logging
from .processor import EidosJsonProcessor, EidosJsonLdProcessor

logger = logging.getLogger('eidos')


try:
    # For text reading
    from .eidos_reader import EidosReader
    eidos_reader = EidosReader()
except Exception as e:
    logger.error('Could not instantiate Eidos reader, text reading '
                 'will not be available.')
    logger.exception(e)
    eidos_reader = None


def process_text(text, out_format='json', save_json='eidos_output.json'):
    """Return an EidosProcessor by processing the given text.

    This constructs a reader object via Java and extracts mentions
    from the text. It then serializes the mentions into JSON and
    processes the result with process_json.

    Parameters
    ----------
    text : str
        The text to be processed.
    out_format : str
        The type of Eidos output to read into and process. Can be one of
        "json" or "json_ld". Default: "json"
    save_json : Optional[str]
        The name of a file in which to dump the JSON output of Eidos.

    Returns
    -------
    ep : EidosJsonProcessor or EidosJsonLdProcessor depending on out_format
        A EidosJsonProcessor or EidosJsonLdProcessor containing the extracted
        INDRA Statements in ep.statements.
    """
    if eidos_reader is None:
        logger.error('Eidos reader is not available.')
        return None
    json_dict = eidos_reader.process_text(text, out_format)
    if save_json:
        with open(save_json, 'wt') as fh:
            json.dump(json_dict, fh, indent=2)
    if out_format == 'json':
        return process_json(json_dict)
    elif out_format == 'json_ld':
        return process_json_ld(json_dict)
    else:
        logger.error('Output format %s is invalid.' % output_format)
        return None


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
    ep : EidosJsonProcessor
        A EidosJsonProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_str(json_str)
    except IOError:
        logger.exception('Could not read file %s.' % file_name)


def process_json_ld_file(file_name):
    """Return an EidosProcessor by processing the given Eidos JSON-LD file.

    The output from the Eidos reader is in json-LD format. This function is
    useful if the output is saved as a file and needs to be processed.

    Parameters
    ----------
    file_name : str
        The name of the JSON-LD file to be processed.

    Returns
    -------
    ep : EidosJsonLdProcessor
        A EidosJsonLdProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_ld_str(json_str)
    except IOError:
        logger.exception('Could not read file %s.' % file_name)


def process_json_str(json_str):
    """Return an EidosProcessor by processing the given Eidos json string.

    The output from the Eidos parser is in json format.

    Parameters
    ----------
    json_str : str
        The json string to be processed.

    Returns
    -------
    ep : EidosJsonProcessor
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


def process_json_ld_str(json_str):
    """Return an EidosJsonLdProcessor by processing the Eidos JSON-LD string.

    The output from the Eidos parser is in JSON-LD format.

    Parameters
    ----------
    json_str : str
        The json-LD string to be processed.

    Returns
    -------
    ep : EidosJsonLdProcessor
        A EidosJsonLdProcessor containing the extracted INDRA Statements
        in ep.statements.
    """
    if not isinstance(json_str, basestring):
        raise TypeError('{} is {} instead of {}'.format(json_str,
                                                        json_str.__class__,
                                                        basestring))
    try:
        json_dict = json.loads(json_str)
    except ValueError:
        logger.error('Could not decode JSON-LD string.')
        return None
    return process_json_ld(json_dict)


def process_json(json_dict):
    """Return an EidosJsonProcessor by processing the given Eidos JSON dict.

    Parameters
    ----------
    json_dict : dict
        The JSON dict to be processed.

    Returns
    -------
    ep : EidosJsonProcessor
        A EidosJsonProcessor containing the extracted INDRA Statements
        in ep.statements.
    """

    ep = EidosJsonProcessor(json_dict)
    ep.get_events()
    return ep


def process_json_ld(json_dict):
    """Return an EidosJsonLdProcessor by processing a Eidos JSON-LD dict.

    Parameters
    ----------
    json_dict : dict
        The JSON-LD dict to be processed.

    Returns
    -------
    ep : EidosJsonLdProcessor
        A EidosJsonLdProcessor containing the extracted INDRA Statements
        in ep.statements.
    """

    ep = EidosJsonLdProcessor(json_dict)
    ep.get_events()
    return ep
