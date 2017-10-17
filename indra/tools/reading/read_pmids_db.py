from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from argparse import ArgumentParser
from docutils.io import InputError
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


def _convert_id_entry(id_entry):
    ret = [s.strip() for s in id_entry.split(':')]
    if len(ret) != 2:
        raise InputError('Improperly formatted id entry: \"%s\"' % id_entry)
    ret[0] = ret[0].lower()
    if ret[0] not in ['pmid', 'doi', 'pmcid']:
        raise InputError('Invalid id type: \"%s\"' % ret[0])
    return ret


def get_content(id_str_list):
    "Load all the content that will be read."
    db = get_primary_db()
    tc_list = db.select_all(db.TextContent, db.TextContent.text_ref_id == db.TextRef.id, getattr(db.TextRef, ))
    return tc_list


def read_content():
    "Apply the readers to the content."


def make_statements():
    "Convert the reader output into statements."


if __name__ == "__main__":
    with open(args.id_file, 'r') as f:
        lines = f.readlines()
    for text_content in get_content(lines):
        reading_result_locs = read_content(text_content)
        