from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import cbio_client

def test_get_cancer_studies():
    study_ids = cbio_client.get_cancer_studies('paad')
    assert(len(study_ids) > 0)
    assert('paad_tcga' in study_ids)

def test_get_cancer_types():
    type_ids = cbio_client.get_cancer_types('lung')
    assert(len(type_ids) > 0)

def test_get_genetic_profiles():
    genetic_profiles = \
        cbio_client.get_genetic_profiles('paad_icgc', 'mutation')
    assert(len(genetic_profiles) > 0)

def test_get_num_sequenced():
    num_case = cbio_client.get_num_sequenced('paad_tcga')
    assert(num_case > 0)

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
                                          ['LOXIMVI_SKIN', 'A101D_SKIN'])
    assert len([x for x in muts]) == 2
    assert 'V600E' in muts['LOXIMVI_SKIN']['BRAF']
    assert 'V600E' in muts['A101D_SKIN']['BRAF']
    assert 'I208V' in muts['LOXIMVI_SKIN']['BRAF']
    assert 'I208V' not in muts['A101D_SKIN']['BRAF']
    assert len(muts['LOXIMVI_SKIN']['AKT1']) == 0
    assert len(muts['A101D_SKIN']['AKT1']) == 0


def test_get_profile_data():
    profile_data = cbio_client.get_profile_data(cbio_client.ccle_study,
                                                ['BRAF', 'PTEN'],
                                                'COPY_NUMBER_ALTERATION',
                                                'all')
    assert profile_data['BT20_BREAST']['PTEN'] == -2
    assert profile_data['BT20_BREAST']['BRAF'] == 1
    assert profile_data['LOXIMVI_SKIN']['PTEN'] == 0
    assert profile_data['LOXIMVI_SKIN']['BRAF'] == 0


def test_get_ccle_cna():
    profile_data = cbio_client.get_ccle_cna(['BRAF'])
    assert profile_data['BT20_BREAST']['BRAF'] == 1
    assert profile_data['LOXIMVI_SKIN']['BRAF'] == 0
