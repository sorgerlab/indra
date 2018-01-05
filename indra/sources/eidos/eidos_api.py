from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str, bytes
from past.builtins import basestring
import json
import logging
from .processor import EidosProcessor

logger = logging.getLogger('eidos')

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
    ep = EidosProcessor(json_dict)
    ep.get_events()
    return ep
