from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from argparse import ArgumentParser
if __name__ == '__main__':
    parser = ArgumentParser(
        description='A tool to read content from the database.'
        )
    parser.add_argument(
        'id_file',
        help='A file containt a list of ids of the form <id_type>:<id>.'
        )
    parser.add_argument(
        '-r', '--readers',
        choice=['reach', 'sparser'],
        help='List of readers to be used.'
        )
    args = parser.parse_args()
from indra.db import get_primary_db


def get_content(id_dict):
    "Load all the content that will be read."


def read_content():
    "Apply the readers to the content."


def make_statements():
    "Convert the reader output into statements."


if __name__ == "__main__":
    with open(args.id_file, 'r') as f:
        lines = f.readlines()
        id_dict = {
            'pmid': []
            }