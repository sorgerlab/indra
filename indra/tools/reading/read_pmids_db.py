"""This module provides essential tools to run reading using indra's own
database. This may also be run as a script; for details run: 
`python read_pmids_db --help`
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import re
import tempfile
import logging
import subprocess
import glob
from argparse import ArgumentParser
from docutils.io import InputError
from datetime import datetime
from indra.sources import reach

logger = logging.getLogger('read_db')
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
    parser.add_argument(
        '-b', '--batch',
        help='Select the number of content entries to be read in each batch.',
        default=1000,
        type=int
        )
    args = parser.parse_args()

from indra.util import unzip_string
from indra.db import get_primary_db, formats
from indra.tools.reading.read_pmids import get_mem_total


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


def get_content(id_str_list, batch_size=1000):
    """Load all the content that will be read."""
    db = get_primary_db()
    clauses = get_clauses(id_str_list, db.TextRef)
    q = db.filter_query(
        db.TextContent,
        db.TextContent.text_ref_id == db.TextRef.id,
        *clauses
        )
    return q.yield_per(batch_size)


class ReachError(Exception):
    pass


class ReachReader(object):
    REACH_MEM = 5  # GB
    MEM_BUFFER = 2  # GB
    name = 'REACH'

    def __init__(self, base_dir, n_proc, force_read, force_fulltext):
        init_msg = 'Running %s with:\n' % self.name
        init_msg += '\n'.join([
            'n_proc=%s' % n_proc,
            'force_read=%s' % force_read,
            'force_fulltext=%s' % force_fulltext
            ])
        logger.info(init_msg)
        self.base_dir = base_dir
        self.exec_path, self.version = self._check_reach_env()
        tmp_dir = tempfile.mkdtemp(
            prefix='read_job_%s' % datetime.now().strftime("%Y%m%d%H%M%S"),
            dir=base_dir
            )
        self.tmp_dir = tmp_dir
        conf_fmt_fname = os.path.join(os.path.dirname(__file__),
                                      'reach_conf_fmt.txt')
        self.conf_file_path = os.path.join(self.tmp_dir, 'indra.conf')
        with open(conf_fmt_fname, 'r') as fmt_file:
            fmt = fmt_file.read()
            with open(self.conf_file_path, 'w') as f:
                f.write(
                    fmt.format(tmp_dir=tmp_dir, num_cores=n_proc)
                    )
        return

    def _check_reach_env(self):
        """Check that the environment supports runnig reach."""
        # Get the path to the reach directory.
        path_to_reach = os.environ.get('REACHPATH', None)
        if path_to_reach is None or not os.path.exists(path_to_reach):
            raise ReachError(
                'Reach path unset or invalid. Check REACHPATH environment var.'
                )
        patt = re.compile('reach-(.*?)\.jar')

        # Find the jar file.
        for fname in os.listdir(path_to_reach):
            m = patt.match(fname)
            if m is not None:
                reach_ex = os.path.join(path_to_reach, fname)
                break
        else:
            raise ReachError("Could not find reach jar in reach dir.")

        logger.info('Using REACH jar at: %s' % reach_ex)

        # Get the reach version.
        reach_version = os.environ.get('REACH_VERSION', None)
        if reach_version is None:
            logger.info('REACH version not set in REACH_VERSION')
            reach_version = re.sub('-SNAP.*?$', '', m.groups()[0])

        logger.info('Using REACH version: %s' % reach_version)
        return reach_ex, reach_version

    def prep_input(self, text_content):
        """Apply the readers to the content."""
        def write_content_file(ext):
            with open('%s.%s' % (text_content.id, ext), 'w') as f:
                f.write(unzip_string(text_content.content))

        if text_content.text_type == formats.XML:
            write_content_file('nxml')
        if text_content.text_type == formats.TEXT:
            write_content_file('txt')
        if text_content.text_type == formats.JSON:
            raise ReachError("I do not know how to handle JSON.")

    def get_output(self):
        """Get the output of a reading job as a list of filenames."""
        json_files = glob.glob(os.path.join(self.tmp_dir, 'output/*.json'))
        json_prefixes = set()
        tc_id_fname_dict = {}
        for json_file in json_files:
            prefix = os.path.basename(json_file).split('.')[0]
            json_prefixes.add(prefix)

        for prefix in json_prefixes:
            tc_id_fname_dict[prefix] = reach.join_json_files(prefix)
        return tc_id_fname_dict

    def convert_output_to_stmts(self):
        pass

    def read(self, tc_list, verbose=False):
        """Read the content.

        This should take in a list of text_content objects, and return nothing.
        """
        mem_tot = get_mem_total()
        if mem_tot is not None and mem_tot <= self.REACH_MEM + self.MEM_BUFFER:
            logger.error(
                "Too little memory to run reach. At least %s required." %
                self.REACH_MEM + self.MEM_BUFFER
                )
            logger.info("REACH not run.")
        elif len(tc_list) > 0:
            # Run REACH!
            logger.info("Beginning reach.")
            args = [
                'java',
                '-Dconfig.file=%s' % self.conf_file_path,
                '-jar', self.exec_path
                ]
            p = subprocess.Popen(args, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            if verbose:
                for line in iter(p.stdout.readline, b''):
                    logger.info(line)
            p_out, p_err = p.communicate()
            if p.returncode:
                logger.error('Problem running REACH:')
                logger.error('Stdout: %s' % p_out.decode('utf-8'))
                logger.error('Stderr: %s' % p_err.decode('utf-8'))
                raise ReachError("Problem running REACH")
        return


def make_statements():
    """Convert the reader output into statements."""


if __name__ == "__main__":
    with open(args.id_file, 'r') as f:
        lines = f.readlines()
    for text_content in get_content(lines):
        pass
