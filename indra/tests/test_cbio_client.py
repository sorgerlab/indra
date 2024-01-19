from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import cbio_client
import pytest

from unittest import skip
@skip('COCI web service not working currently')

@pytest.mark.webservice
def test_get_cancer_studies():
    study_ids = cbio_client.get_cancer_studies('paad')
    assert len(study_ids) > 0
    assert 'paad_tcga' in study_ids


@pytest.mark.webservice
def test_get_cancer_types():
    type_ids = cbio_client.get_cancer_types('lung')
    assert len(type_ids) > 0


@pytest.mark.webservice
def test_get_genetic_profiles():
    genetic_profiles = \
        cbio_client.get_genetic_profiles('paad_icgc', 'mutation')
    assert len(genetic_profiles) > 0


@pytest.mark.webservice
def test_get_num_sequenced():
    num_case = cbio_client.get_num_sequenced('paad_tcga')
    assert num_case > 0


@pytest.mark.webservice
def test_get_ccle_lines_for_mutation():
    """Check how many lines have BRAF V600E mutations.

    Check that this returns a list greater than zero, and more specificially,
    equal to 55 cell lines.
    """
    cl_BRAF_V600E = cbio_client.get_ccle_lines_for_mutation('BRAF', 'V600E')
    assert len(cl_BRAF_V600E) == 55


@pytest.mark.webservice
def test_get_ccle_mutations():
    muts = cbio_client.get_ccle_mutations(['BRAF', 'AKT1'],
                                          ['LOXIMVI_SKIN', 'A101D_SKIN'])
    assert len([x for x in muts]) == 2
    assert 'V600E' in muts['LOXIMVI_SKIN']['BRAF']
    assert 'V600E' in muts['A101D_SKIN']['BRAF']
    assert 'I208V' in muts['LOXIMVI_SKIN']['BRAF']
    assert 'I208V' not in muts['A101D_SKIN']['BRAF']
    assert len(muts['LOXIMVI_SKIN']['AKT1']) == 0
    assert len(muts['A101D_SKIN']['AKT1']) == 0


@pytest.mark.webservice
def test_get_profile_data():
    profile_data = cbio_client.get_profile_data(cbio_client.ccle_study,
                                                ['BRAF', 'PTEN'],
                                                'COPY_NUMBER_ALTERATION',
                                                'all')
    assert profile_data['BT20_BREAST']['PTEN'] == -2
    assert profile_data['BT20_BREAST']['BRAF'] == 1
    assert profile_data['LOXIMVI_SKIN']['PTEN'] == 0
    assert profile_data['LOXIMVI_SKIN']['BRAF'] == 0
    assert len(profile_data) > 0


@pytest.mark.webservice
def test_get_ccle_cna():
    profile_data = cbio_client.get_ccle_cna(['BRAF', 'AKT1'],
                                            ['LOXIMVI_SKIN', 'SKMEL30_SKIN'])
    assert profile_data['SKMEL30_SKIN']['BRAF'] == 1
    assert profile_data['SKMEL30_SKIN']['AKT1'] == -1
    assert profile_data['LOXIMVI_SKIN']['BRAF'] == 0
    assert profile_data['LOXIMVI_SKIN']['AKT1'] == 0
    assert len(profile_data) == 2


@pytest.mark.webservice
def test_get_ccle_mrna():
    mrna = cbio_client.get_ccle_mrna(['XYZ', 'MAP2K1'], ['A375_SKIN'])
    assert 'A375_SKIN' in mrna
    assert mrna['A375_SKIN'] is not None
    assert mrna['A375_SKIN']['MAP2K1'] > 10
    assert mrna['A375_SKIN']['XYZ'] is None
    mrna = cbio_client.get_ccle_mrna(['EGFR', 'BRAF'], ['XXX'])
    assert 'XXX' in mrna
    assert mrna['XXX'] is None


@pytest.mark.webservice
def test_get_ccle_cna_big():
    """
    Get the CNA data on 124 genes in 4 cell lines. Expect to have CNA values
    that are {-2.0, -1.0, 0.0, 1.0, 2.0}. This tests the function at
    a greater scale. Also, test the cell lines' BRAF CNAs
    """
    genes = ["FOSL1", "GRB2", "RPS6KA3", "EIF4EBP1", "DUSP1", "PLXNB1", "SHC2",
             "CBL", "E2F2", "KRAS", "RPS6KA1", "AKT2", "PRKAG2", "JUN", "ELK1",
             "MTOR", "PPP1CA", "TP73", "PIK3R1", "PIK3CG", "FOS", "MLST8",
             "FGFR3", "PRKAG1", "RAF1", "PIK3CA", "TSC1", "AKT3", "SAV1",
             "CCND1", "ETS1", "EXOC7", "CBLB", "MAP2K2", "CDC25A", "PEA15",
             "YAP1", "IRS1", "RPS6KA6", "SOS1", "MDM2", "PIK3R6", "PIK3R2",
             "RASSF1", "AKT1", "BUB1", "PRKAB2", "ETS2", "PRKAG3", "CDK2",
             "TP53", "DAB2IP", "MAPK3", "CDK4", "PRKAB1", "CDK6", "PIK3R5",
             "PAK3", "MYC", "INSR", "SHOC2", "CDKN1A", "STK4", "STK11",
             "RPS6KB1", "RPS6KA2", "KSR1", "PIN1", "MAPK1", "E2F1", "RAC1",
             "CDC6", "EGFR", "NFE2L2", "SHC3", "PTPN11", "PIK3R3", "SHC4",
             "ROCK2", "RPTOR", "CCNA2", "HRAS", "TSC2", "BRCA2", "ALK",
             "MAP2K1", "BRIP1", "CBLC", "RHOA", "RPS6KB2", "FGFR2", "FGFR1",
             "ERBB2", "SOS2", "PRKAA1", "TIAM1", "MET", "PRKAA2", "ROS1",
             "TBK1", "RB1", "PEBP1", "DUSP6", "PTEN", "PDPK1", "MRAS", "NRAS",
             "BRAF", "STK3", "CCNA1", "NFKB1", "CCND3", "PAK4", "KSR2", "ECT2",
             "BRCA1", "VAV1", "CCND2", "ARAF", "PAK2", "PIK3CD", "SHC1",
             "PAK1", "RHEB"]
    cell_lines = ['COLO679_SKIN', 'A2058_SKIN', 'IGR39_SKIN', 'HS294T_SKIN']
    cna = cbio_client.get_ccle_cna(genes, cell_lines)
    values = set()
    for cl in cna:
        for g in cna[cl]:
            values.add(cna[cl][g])
    assert values == {-2.0, -1.0, 0.0, 1.0, 2.0}
    assert cna['COLO679_SKIN']['BRAF'] == 2
    assert cna['A2058_SKIN']['BRAF'] == 1
    assert cna['IGR39_SKIN']['BRAF'] == 1
    assert cna['HS294T_SKIN']['BRAF'] == 1
