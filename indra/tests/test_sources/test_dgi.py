# -*- coding: utf-8 -*-

"""Tests for the DGI processor."""

import os

import pandas as pd

from indra.sources.dgi.processor import DGIProcessor, USECOLS, process_df

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_FILE = os.path.join(HERE, 'resources', 'dgi_sample_interactions.tsv')


def test_dgi_processor():
    df = pd.read_csv(TEST_FILE, sep='\t', usecols=USECOLS, dtype=str)
    df = process_df(df)
    processor = DGIProcessor(df)
    statements = processor.extract_statements()


if __name__ == '__main__':
    test_dgi_processor()
