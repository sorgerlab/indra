# -*- coding: utf-8 -*-

"""Tests for the DGI processor."""

import os

import pandas as pd

from indra.sources.dgi.api import USECOLS, process_df
from indra.sources.dgi.processor import DGIProcessor
from indra.statements import Inhibition

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_FILE = os.path.join(HERE, 'resources', 'dgi_sample_interactions.tsv')


def test_dgi_processor():
    """Test the DGI processor."""
    df = pd.read_csv(TEST_FILE, sep='\t', usecols=USECOLS, dtype=str)
    df = process_df(df)
    processor = DGIProcessor(df)
    statements = processor.extract_statements()
    statement = statements[0]
    assert isinstance(statement, Inhibition)
    assert statement.obj.name == 'CDK7'
    assert statement.obj.db_refs['EGID'] == '1022'
    assert statement.obj.db_refs['HGNC'] == '1778'
    assert statement.subj.name == 'N-(5-\{[(5-tert-butyl-1,3-oxazol-2-yl)methyl]sulfanyl\}-1,3-thiazol-2-yl)piperidine-4-carboxamide', statement.subj.name
    assert statement.subj.db_refs['CHEMBL'] == 'CHEMBL296468'
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {'source': 'CancerCommons', 'interactions': 'inhibitor'}


if __name__ == '__main__':
    test_dgi_processor()
