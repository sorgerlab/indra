import os

from indra.sources.gnbr.processor import *
import indra.sources.gnbr.api as api
from indra.statements.validate import assert_valid_statements


def test_standardize_agent():
    agent = get_std_gene('xxx', '673')
    assert isinstance(agent[0], Agent)
    assert agent[0].name == 'BRAF'
    assert agent[0].db_refs.get('TEXT') == 'xxx'
    assert agent[0].db_refs.get('EGID') == '673'
    assert agent[0].db_refs.get('HGNC') == '1097'


def test_multiple_genes():
    agents = get_std_gene('Erk1/2', '5594;5595')
    assert agents[0].name == 'MAPK1'
    assert agents[1].name == 'MAPK3'
    assert agents[0].db_refs['TEXT'] == 'Erk1/2'
    assert agents[1].db_refs['TEXT'] == 'Erk1/2'
    assert agents[0].db_refs['HGNC'] == '6871'
    assert agents[1].db_refs['HGNC'] == '6877'


def test_process_gene_gene():
    test_path1: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_gene_gene_part1_test.tsv')
    test_path2: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_gene_gene_part2_test.tsv')
    gp = api.process_gene_gene(test_path1, test_path2)
    assert len(gp.statements) != 0
    assert isinstance(gp, GnbrProcessor)
    assert gp.first_type == 'gene'
    assert gp.second_type == 'gene'
    assert isinstance(gp.statements[0], Activation)
    assert isinstance(gp.statements[1], Activation)
    assert isinstance(gp.statements[2], IncreaseAmount)
    assert isinstance(gp.statements[3], IncreaseAmount)
    assert isinstance(gp.statements[4], Complex)
    assert_valid_statements(gp.statements)


def test_process_chemical_gene():
    test_path1: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_chemical_gene_part1_test.tsv')
    test_path2: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_chemical_gene_part2_test.tsv')
    gp = api.process_chemical_gene(test_path1, test_path2)
    assert len(gp.statements) != 0
    assert isinstance(gp, GnbrProcessor)
    assert gp.first_type == 'chemical'
    assert gp.second_type == 'gene'
    assert isinstance(gp.statements[0], Activation)
    assert isinstance(gp.statements[1], Inhibition)
    assert isinstance(gp.statements[2], Complex)
    assert isinstance(gp.statements[3], DecreaseAmount)
    assert_valid_statements(gp.statements)


def test_process_gene_disease():
    test_path1: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_gene_disease_part1_test.tsv')
    test_path2: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_gene_disease_part2_test.tsv')
    gp = api.process_gene_disease(test_path1, test_path2)
    assert len(gp.statements) != 0
    assert isinstance(gp, GnbrProcessor)
    assert gp.first_type == 'gene'
    assert gp.second_type == 'disease'
    assert isinstance(gp.statements[0], Inhibition)
    assert isinstance(gp.statements[1], Activation)
    assert_valid_statements(gp.statements)


def test_process_chemical_disease():
    test_path1: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_chemical_disease_part1_test.tsv')
    test_path2: str = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'gnbr_chemical_disease_part2_test.tsv')
    gp = api.process_chemical_disease(test_path1, test_path2)
    assert len(gp.statements) != 0
    assert isinstance(gp, GnbrProcessor)
    assert gp.first_type == 'chemical'
    assert gp.second_type == 'disease'
    assert isinstance(gp.statements[0], Inhibition)
    assert isinstance(gp.statements[1], Inhibition)
    assert isinstance(gp.statements[2], Inhibition)
    assert isinstance(gp.statements[3], Inhibition)
    assert isinstance(gp.statements[4], Inhibition)
    assert_valid_statements(gp.statements)
