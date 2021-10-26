# -*- coding: utf-8 -*-

"""Tests for the CREEDS processors."""

import pathlib

from indra.sources.creeds import process_from_file
from indra.sources.creeds.processor import LOGGED_MISSING_PART, _get_regulations
from indra.statements import DecreaseAmount, IncreaseAmount, RegulateAmount, BioContext

HERE = pathlib.Path(__file__).parent.resolve()
CREEDS_FOLDER = HERE.joinpath("resources", "creeds_test_data")
GENE_TEST_PATH = CREEDS_FOLDER.joinpath("single_gene.json")
DRUG_TEST_PATH = CREEDS_FOLDER.joinpath("single_drug.json")
DISEASE_TEST_PATH = CREEDS_FOLDER.joinpath("diseases.json")


def test_creeds_gene_processor():
    """Test the CREEDS gene processor."""
    processor = process_from_file(GENE_TEST_PATH, "gene")
    assert 1 == len(processor.records)
    record = processor.records[0]
    assert record["cell_type"] == "heart"
    assert record["geo_id"] == "GSE40601"
    assert record["down_genes"][0][0] == "KRAS"
    assert record["up_genes"][0][0] == "SHC3"
    up_genes, down_genes = _get_regulations(record)
    assert 1 == len(up_genes)
    assert 1 == len(down_genes)

    assert 2 == len(processor.statements), LOGGED_MISSING_PART
    assert all(isinstance(stmt, RegulateAmount) for stmt in processor.statements)
    statement = processor.statements[0]
    assert statement.subj.name == "CDK7"
    assert statement.subj.db_refs["HGNC"] == "1778"
    assert statement.obj.name == "SHC3"
    assert statement.obj.db_refs["HGNC"] == "18181"
    assert isinstance(statement, IncreaseAmount)
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {
        "cell_type": "heart",
        "geo": "GSE40601",
    }
    assert isinstance(evidence.context, BioContext)
    assert evidence.context.species.db_refs["TAXONOMY"] == "9606"

    statement = processor.statements[1]
    assert statement.subj.name == "CDK7"
    assert statement.subj.db_refs["HGNC"] == "1778"
    assert statement.obj.name == "KRAS"
    assert statement.obj.db_refs["HGNC"] == "6407"
    assert isinstance(statement, DecreaseAmount)


def test_creeds_chemical_processor():
    """Test the CREEDS chemical processor."""
    processor = process_from_file(DRUG_TEST_PATH, "chemical")
    assert 1 == len(processor.records)
    record = processor.records[0]
    assert record["cell_type"] == "bone"
    assert record["geo_id"] == "GSE1559"
    assert record["down_genes"][0][0] == "KRAS"
    assert record["up_genes"][0][0] == "SHC3"
    up_genes, down_genes = _get_regulations(record)
    assert 1 == len(up_genes), up_genes
    assert 1 == len(down_genes)
    assert 2 == len(processor.statements), LOGGED_MISSING_PART
    assert all(isinstance(stmt, RegulateAmount) for stmt in processor.statements)
    statement = processor.statements[0]
    assert statement.subj.name == "mitomycin C"
    assert statement.subj.db_refs["DRUGBANK"] == "DB00305"
    assert statement.subj.db_refs["PUBCHEM"] == "5746"
    assert (
        statement.subj.db_refs["SMILES"]
        == "[H][C@]12CN3C4=C([C@@H](COC(N)=O)[C@@]3(OC)[C@@]1([H])N2)C(=O)C(N)=C(C)C4=O"
    )
    assert statement.obj.name == "SHC3"
    assert statement.obj.db_refs["HGNC"] == "18181"
    assert isinstance(statement, IncreaseAmount)
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {
        "cell_type": "bone",
        "geo": "GSE1559",
    }
    assert isinstance(evidence.context, BioContext)
    assert evidence.context.species.db_refs["TAXONOMY"] == "9606"

    statement = processor.statements[1]
    assert statement.subj.name == "mitomycin C"
    assert statement.subj.db_refs["DRUGBANK"] == "DB00305"
    assert statement.subj.db_refs["PUBCHEM"] == "5746"
    assert (
        statement.subj.db_refs["SMILES"]
        == "[H][C@]12CN3C4=C([C@@H](COC(N)=O)[C@@]3(OC)[C@@]1([H])N2)C(=O)C(N)=C(C)C4=O"
    )
    assert statement.obj.name == "KRAS"
    assert statement.obj.db_refs["HGNC"] == "6407"
    assert isinstance(statement, DecreaseAmount)


def test_creeds_disease_processor():
    """Test the CREEDS disease processor."""
    processor = process_from_file(DISEASE_TEST_PATH, "disease")
    assert 1 == len(processor.records)
    record = processor.records[0]
    assert record["cell_type"] == "Muscle - Striated (Skeletal) (MMHCC)"
    assert record["geo_id"] == "GSE466"
    assert record["down_genes"][0][0] == "KRAS"
    assert record["up_genes"][0][0] == "SHC3"
    up_genes, down_genes = _get_regulations(record)
    assert 1 == len(up_genes)
    assert 1 == len(down_genes)

    assert 2 == len(processor.statements), LOGGED_MISSING_PART
    assert all(isinstance(stmt, RegulateAmount) for stmt in processor.statements)
    statement = processor.statements[0]
    assert statement.subj.name == "Shellfish Hypersensitivity"
    assert statement.subj.db_refs["DOID"] == "DOID:0060495"
    assert statement.subj.db_refs["UMLS"] == "C0577625"
    assert statement.obj.name == "SHC3"
    assert statement.obj.db_refs["HGNC"] == "18181"
    assert isinstance(statement, IncreaseAmount)
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {
        "cell_type": "Muscle - Striated (Skeletal) (MMHCC)",
        "geo": "GSE466",
    }
    assert isinstance(evidence.context, BioContext)
    assert evidence.context.species.db_refs["TAXONOMY"] == "9606"

    statement = processor.statements[1]
    assert statement.subj.name == "Shellfish Hypersensitivity"
    assert statement.subj.db_refs["DOID"] == "DOID:0060495"
    assert statement.subj.db_refs["UMLS"] == "C0577625"
    assert statement.obj.name == "KRAS"
    assert statement.obj.db_refs["HGNC"] == "6407"
    assert isinstance(statement, DecreaseAmount)


if __name__ == "__main__":
    test_creeds_gene_processor()
    test_creeds_chemical_processor()
    test_creeds_disease_processor()
