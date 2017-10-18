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
import json
import random
import zlib
from argparse import ArgumentParser
from docutils.io import InputError
from datetime import datetime

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
        choices=['reach', 'sparser'],
        help='List of readers to be used.',
        nargs='+'
        )
    parser.add_argument(
        '-b', '--batch',
        help='Select the number of content entries to be read in each batch.',
        default=1000,
        type=int
        )
    parser.add_argument(
        '-s', '--sample',
        help='To read just a sample of the pmids in the file.',
        type=int
        )
    parser.add_argument(
        '-i', '--in_range',
        help='A range of pmids to be read in the form <start>:<end>.'
        )
    args = parser.parse_args()

from indra.util import unzip_string
from indra.db import get_primary_db, formats, texttypes
from indra.tools.reading.read_pmids import get_mem_total
from indra.sources import reach


def _get_dir(*args):
    dirname = os.path.join(*args)
    if os.path.isabs(dirname):
        dirpath = dirname
    elif os.path.exists(dirname):
        dirpath = os.path.abspath(dirname)
    else:
        dirpath = os.path.join(os.path.dirname(__file__), dirname)
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    return dirpath


def _time_stamp():
    return datetime.now().strftime("%Y%m%d%H%M%S")


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
    logger.info('Found %d content entries.' % q.count())
    for _, texttype in texttypes.iterattrs():
        logger.info(
            '%d were %ss' % (
                q.filter(db.TextContent.text_type == texttype).count(),
                texttype
                )
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
        self.base_dir = _get_dir(base_dir)
        self.exec_path, self.version = self._check_reach_env()
        tmp_dir = tempfile.mkdtemp(
            prefix='read_job_%s' % _time_stamp(),
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
        self.input_dir = _get_dir(tmp_dir, 'input')
        self.output_dir = _get_dir(tmp_dir, 'output')
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

    def write_content(self, text_content):
        def write_content_file(ext):
            fname = '%s.%s' % (text_content.id, ext)
            with open(os.path.join(self.input_dir, fname), 'w') as f:
                f.write(
                    zlib.decompress(
                        text_content.content, 16+zlib.MAX_WBITS
                        ).decode('utf8')
                    )
            logger.info('%s saved for reading by reach.' % fname)
        if text_content.format == formats.XML:
            write_content_file('nxml')
        elif text_content.format == formats.TEXT:
            write_content_file('txt')
        elif text_content.format == formats.JSON:
            raise ReachError("I do not know how to handle JSON.")
        else:
            raise ReachError("Unrecognized format %s." % text_content.format)

    def prep_input(self, tc_list):
        """Apply the readers to the content."""
        logger.info("Prepping input.")
        for text_content in tc_list:
            self.write_content(text_content)
        return

    def get_output(self):
        """Get the output of a reading job as a list of filenames."""
        logger.info("Getting outputs.")
        # Get the set of prefixes (each will correspond to three json files.)
        json_files = glob.glob(os.path.join(self.output_dir, '*.json'))
        json_prefixes = set()
        for json_file in json_files:
            prefix = os.path.basename(json_file).split('.')[0]
            json_prefixes.add(os.path.join(self.output_dir, prefix))

        # Join each set of json files and store the json dict.
        tc_id_fname_dict = {}
        for prefix in json_prefixes:
            tc_id_fname_dict[prefix] = reach.join_json_files(prefix)
            logger.info('Joined files for prefix %s.' % prefix)
        return tc_id_fname_dict

    def convert_output_to_stmts(self, json_dict):
        """Process the REACH output with INDRA"""
        # Convert the JSON object into a string first so that a series of
        # string replacements can happen in the REACH processor
        reach_proc = reach.process_json_str(json.dumps(json_dict))
        return reach_proc.statements

    def read(self, tc_list, verbose=False):
        """Read the content.

        This should take in a list of text_content objects, and return nothing.
        """
        ret = None
        mem_tot = get_mem_total()
        if mem_tot is not None and mem_tot <= self.REACH_MEM + self.MEM_BUFFER:
            logger.error(
                "Too little memory to run reach. At least %s required." %
                self.REACH_MEM + self.MEM_BUFFER
                )
            logger.info("REACH not run.")
        elif len(tc_list) > 0:
            # Prep the content
            self.prep_input(tc_list)
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
            logger.info("Reach finished.")
            ret = self.get_output()
        return ret


def read(tc_list, readers, *args, **kwargs):
    """Perform the reading, returning dicts of jsons."""
    output_dict = {}
    if 'reach' in readers:
        r = ReachReader(*args, **kwargs)
        output_dict.update(r.read(tc_list))
        logger.info("Read %d text content entries with reach."
                    % len(output_dict))
    return output_dict


def post_reading_output():
    """Put the reading output on the database."""


def make_statements():
    """Convert the reader output into statements."""


if __name__ == "__main__":
    base_dir = _get_dir('run_%s' % ('_and_'.join(args.readers)))

    # Process the arguments.
    with open(args.id_file, 'r') as f:
        lines = f.readlines()

    if args.sample is not None:
        lines = random.sample(lines, args.sample)

    if args.in_range is not None:
        start_idx, end_idx = [int(n) for n in args.in_range.split(':')]
        lines = lines[start_idx:end_idx]

    # Start reading in batches.
    batch_list = []
    for i, text_content in enumerate(get_content(lines, args.batch)):
        batch_list.append(text_content)
        if (i+1) % args.batch is 0:
            read(batch_list, base_dir, 1, True, False)
            batch_list = []
    if len(batch_list) > 0:
        read(batch_list, base_dir, 1, True, False)
