from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import ndex2.client
from .processor import NdexCxProcessor


def process_cx_file(file_name, require_grounding=True):
    """Process a CX JSON file into Statements.

    Parameters
    ----------
    file_name : str
        Path to file containing CX JSON.
    require_grounding: bool
        Whether network nodes lacking grounding information should be included
        among the extracted Statements (default is True).

    Returns
    -------
    NdexCxProcessor
        Processor containing Statements.
    """
    with open(file_name, 'rt') as fh:
        json_list = json.load(fh)
        return process_cx(json_list, require_grounding=require_grounding)


def process_ndex_network(network_id, username=None, password=None,
                         require_grounding=True):
    """Process an NDEx network into Statements.

    Parameters
    ----------
    network_id : str
        NDEx network ID.
    username : str
        NDEx username.
    password : str
        NDEx password.
    require_grounding: bool
        Whether network nodes lacking grounding information should be included
        among the extracted Statements (default is True).

    Returns
    -------
    NdexCxProcessor
        Processor containing Statements. Returns None if there if the HTTP
        status code indicates an unsuccessful request.
    """
    nd = ndex2.client.Ndex2(username=username, password=password)
    res = nd.get_network_as_cx_stream(network_id)
    if res.status_code != 200:
        logger.error('Problem downloading network: status code %s' %
                     res.status_code)
        logger.error('Response: %s' % res.text)
        return None
    json_list = res.json()
    return process_cx(json_list, require_grounding=require_grounding)


def process_cx(cx_json, require_grounding=True):
    """Process a CX JSON object into Statements.

    Parameters
    ----------
    cx_json : list
        CX JSON object.
    require_grounding: bool
        Whether network nodes lacking grounding information should be included
        among the extracted Statements (default is True).

    Returns
    -------
    NdexCxProcessor
        Processor containing Statements.
    """
    ncp = NdexCxProcessor(cx_json, require_grounding=require_grounding)
    ncp.get_statements()
    return ncp
