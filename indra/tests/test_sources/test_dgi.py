# -*- coding: utf-8 -*-

"""Tests for the DGI processor."""

import os

import pandas as pd

from indra.sources.dgi.api import USECOLS
from indra.sources.dgi import process_df
from indra.statements import Inhibition

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_FILE = os.path.join(HERE, 'resources', 'dgi_sample_interactions.tsv')


def test_dgi_processor():
    """Test the DGI processor."""
    df = pd.read_csv(TEST_FILE, sep='\t', usecols=USECOLS, dtype=str)
    df = df[USECOLS]
    dp = process_df(df)
    statement = dp.statements[0]
    assert isinstance(statement, Inhibition)
    assert statement.obj.name == 'CDK7'
    assert statement.obj.db_refs['EGID'] == '1022'
    assert statement.obj.db_refs['HGNC'] == '1778'
    assert statement.subj.name == 'BMS-387032'
    assert statement.subj.db_refs['IUPHAR.LIGAND'] == '5670'
    assert 1 == len(statement.evidence)
    evidence = statement.evidence[0]
    assert evidence.pmid is None
    assert evidence.annotations == {'source': 'ChEMBL',
                                    'interactions': 'inhibitor'}


if __name__ == '__main__':
    test_dgi_processor()
