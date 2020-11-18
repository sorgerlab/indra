import json
import time
import openpyxl
import requests
from indra.config import get_config
from .processor import SofiaJsonProcessor, SofiaExcelProcessor


def process_table(fname, extract_filter=None):
    """Return processor by processing a given sheet of a spreadsheet file.

    Parameters
    ----------
    fname : str
        The name of the Excel file (typically .xlsx extension) to process
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute.
    """
    book = openpyxl.load_workbook(fname, read_only=True)
    try:
        rel_sheet = book['Relations']
    except Exception as e:
        rel_sheet = book['Causal']
    event_sheet = book['Events']
    entities_sheet = book['Entities']

    sp = SofiaExcelProcessor(rel_sheet.rows, event_sheet.rows,
                             entities_sheet.rows)
    if extract_filter is None or 'influence' in extract_filter:
        sp.extract_relations(rel_sheet.rows)
    if extract_filter is None or 'event' in extract_filter:
        sp.extract_events(event_sheet.rows, rel_sheet.rows)
    return sp


def process_text(text, out_file='sofia_output.json', auth=None,
                 extract_filter=None):
    """Return processor by processing text given as a string.

    Parameters
    ----------
    text : str
        A string containing the text to be processed with Sofia.
    out_file : Optional[str]
        The path to a file to save the reader's output into.
        Default: sofia_output.json
    auth : Optional[list]
        A username/password pair for the Sofia web service. If not given,
        the SOFIA_USERNAME and SOFIA_PASSWORD values are loaded from either
        the INDRA config or the environment.
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute. If the API did not process
        the text, None is returned.
    """
    text_json = {'text': text}
    if not auth:
        user, password = _get_sofia_auth()
    else:
        user, password = auth
    if not user or not password:
        raise ValueError('Could not use SOFIA web service since'
                         ' authentication information is missing. Please'
                         ' set SOFIA_USERNAME and SOFIA_PASSWORD in the'
                         ' INDRA configuration file or as environmental'
                         ' variables.')
    json_response, status_code, process_status = \
        _text_processing(text_json=text_json, user=user, password=password)

    # Check response status
    if process_status != 'Done' or status_code != 200:
        return None

    # Cache reading output
    if out_file:
        with open(out_file, 'w') as fh:
            json.dump(json_response, fh, indent=1)

    return process_json(json_response, extract_filter=extract_filter)


def process_json(json_obj, extract_filter=None):
    """Return processor by processing a JSON object returned by Sofia.

    Parameters
    ----------
    json_obj : json
        A JSON object containing extractions from Sofia.
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute.
    """
    sp = SofiaJsonProcessor(json_obj)
    if extract_filter is None or 'influence' in extract_filter:
        sp.extract_relations(json_obj)
    if extract_filter is None or 'event' in extract_filter:
        sp.extract_events(json_obj)
    return sp


def process_json_file(fname, extract_filter=None):
    """Return processor by processing a JSON file produced by Sofia.

    Parameters
    ----------
    fname : str
        The name of the JSON file to process
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence' and 'event'. If not given, all relation
        types are extracted. This argument can be used if, for instance,
        only Influence statements are of interest. Default: None

    Returns
    -------
    indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute.
    """
    with open(fname, 'r') as fh:
        jd = json.load(fh)
    return process_json(jd, extract_filter=extract_filter)


def _get_sofia_auth():
    sofia_username = get_config('SOFIA_USERNAME')
    sofia_password = get_config('SOFIA_PASSWORD')
    return sofia_username, sofia_password


def _sofia_api_post(api, option, json, auth):
    return requests.post(url=api + option, json=json, auth=auth)


def _text_processing(text_json, user, password):
    assert len(text_json) > 0

    sofia_api = 'https://sofia.worldmodelers.com'
    auth = (user, password)

    # Initialize process
    resp = _sofia_api_post(api=sofia_api, option='/process_text',
                           json=text_json, auth=auth)
    res_json = resp.json()

    # Get status
    status = _sofia_api_post(api=sofia_api, option='/status',
                             json=res_json, auth=auth)

    # Check status every two seconds
    while status.json()['Status'] == 'Processing':
        time.sleep(2.0)
        status = _sofia_api_post(api=sofia_api, option='/status',
                                 json=res_json, auth=auth)

    results = _sofia_api_post(api=sofia_api, option='/results',
                              json=res_json, auth=auth)
    status_code = results.status_code
    process_status = status.json()['Status']
    return results.json(), status_code, process_status
