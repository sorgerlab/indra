# -*- coding: utf-8 -*-

"""
Module that contains the command line app

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will cause
problems--the code will get executed twice:
 - When you run `python3 -m indra` python will execute
   ``__main__.py`` as a script. That means there won't be any
   ``indra.__main__`` in ``sys.modules``.
 - When you import __main__ it will get executed again (as a module) because
   there's no ``indra.__main__`` in ``sys.modules``.
Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import sys

import click
import os

from indra.machine.utils import copy_default_config


@click.group()
def main():
    """INDRA"""


@main.group()
def machine():
    """RAS Machine"""


@machine.command()
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
    os.mkdir(os.path.join(directory, 'json'))
    os.mkdir(os.path.join(directory, 'json', 'abstract'))
    os.mkdir(os.path.join(directory, 'json', 'full'))
    copy_default_config(os.path.join(directory, 'config.yaml'))


if __name__ == '__main__':
    main()
