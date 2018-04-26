"""
Module that contains the command line app

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will cause
problems--the code will get executed twice:
 - When you run `python3 -m indra.tools.machine` python will execute
   ``__main__.py`` as a script. That means there won't be any
   ``indra.tools.machine__main__`` in ``sys.modules``.
 - When you import __main__ it will get executed again (as a module) because
   there's no ``indra.tools.machine__main__`` in ``sys.modules``.
Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import os
import sys

import click

from indra.tools.machine.config import copy_default_config


@click.group()
def main():
    """INDRA"""


@main.command()
@click.argument('directory')
def make(directory):
    """Makes a RAS Machine directory"""

    if os.path.exists(directory):
        if os.path.isdir(directory):
            click.echo('Directory already exists')
        else:
            click.echo('Path exists and is not a directory')
        sys.exit()

    os.makedirs(directory)
    os.mkdir(os.path.join(directory, 'jsons'))
    copy_default_config(os.path.join(directory, 'config.yaml'))


@main.command()
@click.argument('model_path')
@click.option('--config', help='Specify configuration file path, otherwise '
                               'looks for config.yaml in model path')
@click.option('-d', '--num-days', type=int)
def run_with_search(model_path, config, num_days):
    """Run with PubMed search for new papers."""
    from indra.tools.machine.machine import run_with_search_helper
    run_with_search_helper(model_path, config, num_days=num_days)


@main.command()
@click.argument('model_path')
def summarize(model_path):
    """Print model summary."""
    from indra.tools.machine.machine import summarize_helper
    summarize_helper(model_path)


@main.command()
@click.argument('model_path')
@click.option('--pmids', type=click.File(), default=sys.stdin,
              help="A file with a PMID on each line")
def run_with_pmids(model_path, pmids):
    """Run with given list of PMIDs."""
    from indra.tools.machine.machine import run_with_pmids_helper
    run_with_pmids_helper(model_path, pmids)


@main.command()
@click.argument('model_path')
@click.argument('name')
@click.option('--output', type=click.File('w'))
def to_cx(model_path, name, output):
    from indra.tools.machine.machine import load_model, assemble_cx
    model = load_model(model_path)
    stmts = model.get_statements()
    cx_str = assemble_cx(stmts, name)
    click.echo(cx_str, output)


if __name__ == '__main__':
    main()
