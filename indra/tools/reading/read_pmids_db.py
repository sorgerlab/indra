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


def _convert_id_entry(id_entry, allowed_types=None):
    ret = [s.strip() for s in id_entry.split(':')]
    if len(ret) != 2:
        raise InputError('Improperly formatted id entry: \"%s\"' % id_entry)
    ret[0] = ret[0].lower()
    if allowed_types is not None and ret[0] not in allowed_types:
        raise InputError('Invalid id type: \"%s\"' % ret[0])
    return ret


def get_clauses(id_str_list, table):
    """Get a list of clauses to be passed to a db query."""
    id_types = table.__table__.columns.keys()
    id_dict = {id_type: [] for id_type in id_types}
    for id_entry in id_str_list:
        id_type, id_val = _convert_id_entry(id_entry, id_types)
        id_dict[id_type].append(id_val)
    return [getattr(table, id_type).in_(id_list)
            for id_type, id_list in id_dict.items() if len(id_list)]


def get_content(id_str_list, batch_size = 1000):
    """Load all the content that will be read."""
    db = get_primary_db()
    clauses = get_clauses(id_str_list, db.TextRef)
    q = db.filter_query(
        db.TextContent,
        db.TextContent.text_ref_id == db.TextRef.id,
        *clauses
        )
    return q.yield_per(batch_size)


def read_content():
    """Apply the readers to the content."""


def make_statements():
    """Convert the reader output into statements."""


if __name__ == "__main__":
    with open(args.id_file, 'r') as f:
        lines = f.readlines()
    for text_content in get_content(lines):
        reading_result_locs = read_content(text_content)
