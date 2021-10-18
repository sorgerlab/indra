# -*- coding: utf-8 -*-

"""Tests for the CREEDS processors."""

import pathlib

from indra.sources.creeds import process_from_file
from indra.statements import DecreaseAmount, RegulateAmount

HERE = pathlib.Path(__file__).parent.resolve()
CREEDS_FOLDER = HERE.joinpath("resources", "creeds_test_data")
GENE_TEST_PATH = CREEDS_FOLDER.joinpath("single_gene.json")


def test_creeds_processor():
    """Test the CREEDS processor."""
    processor = process_from_file(GENE_TEST_PATH, "gene")
    assert 4 == len(processor.statements)
    assert all(isinstance(stmt, RegulateAmount) for stmt in processor)

    statement = processor.statements[0]
    assert isinstance(statement, DecreaseAmount)
    assert statement.obj.name == 'Pdk4'


if __name__ == '__main__':
    test_creeds_processor()
