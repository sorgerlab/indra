# -*- coding: utf-8 -*-

"""Tests for the CREEDS processors."""

import pathlib

from indra.sources.creeds import process_from_file
from indra.sources.creeds.processor import LOGGED_MISSING_PART, _get_regulations
from indra.statements import DecreaseAmount, RegulateAmount

HERE = pathlib.Path(__file__).parent.resolve()
CREEDS_FOLDER = HERE.joinpath("resources", "creeds_test_data")
GENE_TEST_PATH = CREEDS_FOLDER.joinpath("single_gene.json")


def test_creeds_processor():
    """Test the CREEDS processor."""
    processor = process_from_file(GENE_TEST_PATH, "gene")
    assert 1 == len(processor.records)
    record = processor.records[0]
    assert record["cell_type"] == "heart"
    assert record["geo_id"] == "GSE44192"
    assert record["mm_gene_symbol"] == "Plin5"
    assert record["down_genes"][0][0] == "Pdk4"
    up_genes, down_genes = _get_regulations(record)
    assert 2 == len(up_genes)
    assert 2 == len(down_genes)

    assert 4 == len(processor.statements), LOGGED_MISSING_PART
    assert all(isinstance(stmt, RegulateAmount) for stmt in processor.statements)
    statement = processor.statements[0]
    assert isinstance(statement, DecreaseAmount)
    assert statement.subj.name == "Plin5"
    assert statement.subj.db_refs["MGI"] == "1914218"
    assert statement.obj.name == "Pdk4"
    assert statement.obj.db_refs["MGI"] == "1351481"
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {
        "organism": "mouse",
        "cell": "heart",
        "geo": "GSE44192",
    }


if __name__ == "__main__":
    test_creeds_processor()
