import time

import openpyxl
import requests

from .processor import SofiaProcessor


def process_table(fname):
    """Return processor by processing a given sheet of a spreadsheet file.

    Parameters
    ----------
    fname : str
        The name of the Excel file (typically .xlsx extension) to process

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute
    """
    book = openpyxl.load_workbook(fname, read_only=True)
    try:
        rel_sheet = book['Relations']
    except Exception as e:
        rel_sheet = book['Causal']
    event_sheet = book['Events']
    entities_sheet = book['Entities']
    sp = SofiaProcessor(rel_sheet.rows, event_sheet.rows, entities_sheet.rows)
    return sp


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
        if status.json()['Status'] == 'Done':
            # Get results when processing is done
            results = _sofia_api_post(api=sofia_api, option='/results',
                                      json=res_json, auth=auth)
            return results.json()

    # The while loop exited without 'Done' status; return the api response
    results = _sofia_api_post(api=sofia_api, option='/results',
                              json=res_json, auth=auth)
    return results.json()
