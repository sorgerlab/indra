from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import cbio_client


def test_send_request_ccle():
    """Sends a request and gets back a dataframe of all cases in ccle study.

    Check that the dataframe is longer than one.
    """
    data = {'cmd': 'getCaseLists',
            'cancer_study_id': 'cellline_ccle_broad'}
    df = cbio_client.send_request(**data)
    assert(len(df) > 0)


def test_get_ccle_lines_for_mutation():
    """Check how many lines have BRAF V600E mutations.

    Check that this returns a list greater than zero, and more specificially,
    equal to 55 cell lines.
    """
    cl_BRAF_V600E = cbio_client.get_ccle_lines_for_mutation('BRAF', 'V600E')
    assert(len(cl_BRAF_V600E) == 55)


def test_get_mutations_ccle():
    muts = cbio_client.get_mutations_ccle(['BRAF', 'AKT1'],
                                          ['LOXIMVI', 'A101D'])
    assert len([x for x in muts]) == 2
    assert 'V600E' in muts['LOXIMVI']['BRAF']
    assert 'V600E' in muts['A101D']['BRAF']
    assert 'I208V' in muts['LOXIMVI']['BRAF']
    assert 'I208V' not in muts['A101D']['BRAF']
    assert len(muts['LOXIMVI']['AKT1']) == 0
    assert len(muts['A101D']['AKT1']) == 0
