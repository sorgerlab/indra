# -*- coding: utf-8 -*-

import sys

import click

from indra.machine.utils import run_with_pmids_helper, run_with_search_helper, summarize_helper


@click.group()
def main():
    """The RAS Machine and utilities"""


@main.command()
@click.argument('model_path')
@click.option('--config', help='Specify configuration file path, otherwise '
                               'looks for config.yaml in model path')
def run_with_search(model_path, config):
    """Run with PubMed search for new papers."""
    run_with_search_helper(model_path, config)


@main.command()
@click.argument('model_path')
def summarize(model_path):
    """Print model summary."""
    summarize_helper(model_path)


@main.command()
@click.argument('model_path')
@click.option('--pmids', type=click.File(), default=sys.stdin,
              help="A file with a PMID on each line")
def run_with_pmids(model_path, pmids):
    """Run with given list of PMIDs."""
    run_with_pmids_helper(model_path, pmids)


if __name__ == '__main__':
    main()
