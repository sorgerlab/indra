from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import cbio_client


def test_send_request_ccle():
    '''
    Sends a request and gets back a dataframe of all cases in ccle study
    Check that the dataframe is longer than one
    '''
    data = {'cmd': 'getCaseLists',
        'cancer_study_id': 'cellline_ccle_broad'}
    df = cbio_client.send_request(data)
    assert(len(df) > 0)
