# -*- coding: utf-8 -*-

import sys

import click


@click.group()
def machine():
    """The RAS Machine and utilities"""


@machine.command()
@click.argument('model_path')
@click.option('--config', help='Specify configuration file path, otherwise '
                               'looks for config.yaml in model path')
def run_with_search(model_path, config):
    """Run with PubMed search for new papers."""
    from indra.tools.machine.utils import run_with_search_helper
    run_with_search_helper(model_path, config)


@machine.command()
@click.argument('model_path')
def summarize(model_path):
    """Print model summary."""
    from indra.tools.machine.utils import summarize_helper
    summarize_helper(model_path)


@machine.command()
@click.argument('model_path')
@click.option('--pmids', type=click.File(), default=sys.stdin,
              help="A file with a PMID on each line")
def run_with_pmids(model_path, pmids):
    """Run with given list of PMIDs."""
    from indra.tools.machine.utils import run_with_pmids_helper
    run_with_pmids_helper(model_path, pmids)


if __name__ == '__main__':
    machine()
