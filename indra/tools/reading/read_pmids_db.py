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
import pickle
import shutil
from argparse import ArgumentParser
from docutils.io import InputError
from datetime import datetime
from math import log10, floor
from os.path import join as pjoin
from os import path
import itertools

logger = logging.getLogger('read_db')
if __name__ == '__main__':
    parser = ArgumentParser(
        description='A tool to read content from the database.'
        )
    parser.add_argument(
        '-i', '--id_file',
        help='A file containt a list of ids of the form <id_type>:<id>.'
        )
    parser.add_argument(
        '-f', '--file_file',
        help=('A file containing a list of files to be input into reach. These'
              'should be nxml or txt files.')
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
        '-I', '--in_range',
        help='A range of pmids to be read in the form <start>:<end>.'
        )
    parser.add_argument(
        '-v', '--verbose',
        help='Include output from the readers.',
        action='store_true'
        )
    parser.add_argument(
        '-q', '--quiet',
        help='Suppress most output. Overrides -v and -d options.',
        action='store_true'
        )
    parser.add_argument(
        '-d', '--debug',
        help='Set the logging to debug level.',
        action='store_true'
        )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

from indra.util import unzip_string
from indra.db import get_primary_db, formats, texttypes
from indra.tools.reading.read_pmids import get_mem_total
from indra.sources import reach


def _get_dir(*args):
    dirname = pjoin(*args)
    if path.isabs(dirname):
        dirpath = dirname
    elif path.exists(dirname):
        dirpath = path.abspath(dirname)
    else:
        dirpath = pjoin(path.dirname(__file__), dirname)
    if not path.exists(dirpath):
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


def get_table_string(q, db):
    """Create a table with some summary data for a query."""
    N_tot = q.count()
    log_n = floor(log10(N_tot))
    cols = list(formats.values()) + ['tot']
    col_fmt = ' %%%ds' % max(4, log_n)
    cols_strs = [col_fmt % fmt for fmt in cols]
    ret_str = 'Summary Statisitics:\n' + ' '*10 + ''.join(cols_strs) + '\n'
    col_counts = dict.fromkeys(formats.values())
    col_counts['tot'] = []
    for texttype in texttypes.values():
        line = '%8s: ' % texttype
        counts = []
        for text_format in cols[:-1]:
            if col_counts[text_format] is None:
                col_counts[text_format] = []
            c = q.filter(
                db.TextContent.text_type == texttype,
                db.TextContent.format == text_format
                ).count()
            line += col_fmt % c
            counts.append(c)
            col_counts[text_format].append(c)
        line += col_fmt % sum(counts)
        ret_str += line + '\n'
        col_counts['tot'].append(sum(counts))
    ret_str += '%8s: ' % 'total' + ''.join([col_fmt % sum(col_counts[col])
                                            for col in cols])
    return ret_str


def get_content(id_str_list, batch_size=1000):
    """Load all the content that will be read."""
    db = get_primary_db()
    logger.debug("Got db handle.")
    clauses = get_clauses(id_str_list, db.TextRef)
    logger.debug("Generated %d clauses." % len(clauses))
    if len(clauses):
        q = db.filter_query(
            db.TextContent,
            db.TextContent.text_ref_id == db.TextRef.id,
            *clauses
            )
        logger.info(get_table_string(q, db))
        ret = q.yield_per(batch_size)
    else:
        logger.info("Did not retreive content from database.")
        ret = []
    return ret


class ReachError(Exception):
    pass


class ReachReader(object):
    REACH_MEM = 5  # GB
    MEM_BUFFER = 2  # GB
    name = 'REACH'

    def __init__(self, base_dir=None, n_proc=1):
        if base_dir is None:
            base_dir = 'run_reach'
        self.n_proc = n_proc
        self.base_dir = _get_dir(base_dir)
        self.exec_path, self.version = self._check_reach_env()
        tmp_dir = tempfile.mkdtemp(
            prefix='read_job_%s' % _time_stamp(),
            dir=base_dir
            )
        self.tmp_dir = tmp_dir
        conf_fmt_fname = pjoin(path.dirname(__file__),
                               'reach_conf_fmt.txt')
        self.conf_file_path = pjoin(self.tmp_dir, 'indra.conf')
        with open(conf_fmt_fname, 'r') as fmt_file:
            fmt = fmt_file.read()
            loglevel = 'INFO'  # 'DEBUG' if logger.level == logging.DEBUG else 'INFO'
            with open(self.conf_file_path, 'w') as f:
                f.write(
                    fmt.format(tmp_dir=tmp_dir, num_cores=n_proc,
                               loglevel=loglevel)
                    )
        self.input_dir = _get_dir(tmp_dir, 'input')
        self.output_dir = _get_dir(tmp_dir, 'output')
        return

    def _check_reach_env(self):
        """Check that the environment supports runnig reach."""
        # Get the path to the reach directory.
        path_to_reach = os.environ.get('REACHPATH', None)
        if path_to_reach is None or not path.exists(path_to_reach):
            raise ReachError(
                'Reach path unset or invalid. Check REACHPATH environment var.'
                )
        patt = re.compile('reach-(.*?)\.jar')

        # Find the jar file.
        for fname in os.listdir(path_to_reach):
            m = patt.match(fname)
            if m is not None:
                reach_ex = pjoin(path_to_reach, fname)
                break
        else:
            raise ReachError("Could not find reach jar in reach dir.")

        logger.debug('Using REACH jar at: %s' % reach_ex)

        # Get the reach version.
        reach_version = os.environ.get('REACH_VERSION', None)
        if reach_version is None:
            logger.debug('REACH version not set in REACH_VERSION')
            reach_version = re.sub('-SNAP.*?$', '', m.groups()[0])

        logger.debug('Using REACH version: %s' % reach_version)
        return reach_ex, reach_version

    def write_content(self, text_content):
        def write_content_file(ext):
            fname = '%s.%s' % (text_content.id, ext)
            with open(pjoin(self.input_dir, fname), 'wb') as f:
                f.write(
                    zlib.decompress(
                        text_content.content, 16+zlib.MAX_WBITS
                        )
                    )
            logger.debug('%s saved for reading by reach.' % fname)
        if text_content.format == formats.XML:
            write_content_file('nxml')
        elif text_content.format == formats.TEXT:
            write_content_file('txt')
        elif text_content.format == formats.JSON:
            raise ReachError("I do not know how to handle JSON.")
        else:
            raise ReachError("Unrecognized format %s." % text_content.format)

    def prep_input(self, read_list):
        """Apply the readers to the content."""
        logger.info("Prepping input.")
        for text_content in read_list:
            if isinstance(text_content, str):
                fname = text_content.strip()
                shutil.copy(
                    fname,
                    pjoin(self.input_dir, path.basename(fname))
                    )
            else:
                self.write_content(text_content)
        return

    def get_output(self):
        """Get the output of a reading job as a list of filenames."""
        logger.info("Getting outputs.")
        # Get the set of prefixes (each will correspond to three json files.)
        json_files = glob.glob(pjoin(self.output_dir, '*.json'))
        json_prefixes = set()
        for json_file in json_files:
            prefix = path.basename(json_file).split('.')[0]
            json_prefixes.add(pjoin(self.output_dir, prefix))

        # Join each set of json files and store the json dict.
        tc_id_fname_dict = {}
        for prefix in json_prefixes:
            base_prefix = path.basename(prefix)
            tc_id_fname_dict[base_prefix] = reach.join_json_files(prefix)
            logger.debug('Joined files for prefix %s.' % base_prefix)
        return tc_id_fname_dict

    def convert_output_to_stmts(self, json_dict):
        """Process the REACH output with INDRA"""
        # Convert the JSON object into a string first so that a series of
        # string replacements can happen in the REACH processor
        reach_proc = reach.process_json_str(json.dumps(json_dict))
        return reach_proc.statements

    def read(self, read_list, verbose=False, force_read=True,
             force_fulltext=False):
        """Read the content."""
        init_msg = 'Running %s with:\n' % self.name
        init_msg += '\n'.join([
            'n_proc=%s' % self.n_proc,
            'force_read=%s' % force_read,
            'force_fulltext=%s' % force_fulltext
            ])
        logger.info(init_msg)
        ret = None
        mem_tot = get_mem_total()
        if mem_tot is not None and mem_tot <= self.REACH_MEM + self.MEM_BUFFER:
            logger.error(
                "Too little memory to run reach. At least %s required." %
                self.REACH_MEM + self.MEM_BUFFER
                )
            logger.info("REACH not run.")
        elif len(read_list) > 0:
            # Prep the content
            self.prep_input(read_list)
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
                    logger.info('REACH: ' + line.strip().decode('utf8'))
            p_out, p_err = p.communicate()
            if p.returncode:
                logger.error('Problem running REACH:')
                logger.error('Stdout: %s' % p_out.decode('utf-8'))
                logger.error('Stderr: %s' % p_err.decode('utf-8'))
                raise ReachError("Problem running REACH")
            logger.info("Reach finished.")
            ret = self.get_output()
        return ret


def read(read_list, readers, *args, **kwargs):
    """Perform the reading, returning dicts of jsons."""
    base_dir = _get_dir('run_%s' % ('_and_'.join(readers)))
    output_dict = {}
    if 'reach' in readers:
        n_proc = kwargs.pop('n_proc', 1)
        r = ReachReader(n_proc=n_proc, base_dir=base_dir)
        output_dict.update(r.read(read_list, *args, **kwargs))
        logger.info("Read %d text content entries with reach."
                    % len(output_dict))
    return output_dict


def post_reading_output():
    """Put the reading output on the database."""


def make_statements():
    """Convert the reader output into statements."""


if __name__ == "__main__":
    # Process the arguments.
    id_lines = []
    if args.id_file is not None:
        with open(args.id_file, 'r') as f:
            id_lines = f.readlines()

    file_lines = []
    if args.file_file is not None:
        with open(args.file_file, 'r') as f:
            file_lines = f.readlines()
    assert id_lines or file_lines, 'No inputs provided.'

    if args.sample is not None:
        id_lines = random.sample(id_lines, args.sample)

    if args.in_range is not None:
        start_idx, end_idx = [int(n) for n in args.in_range.split(':')]
        id_lines = id_lines[start_idx:end_idx]

    # Start reading in batches.
    batch_list = []
    outputs = {}
    read_args = [args.readers]
    read_kwargs = dict(
        verbose=args.verbose and not args.quiet,
        force_read=True,
        force_fulltext=False,
        n_proc=1
        )
    logger.debug("Getting iterator.")
    input_iter = itertools.chain(file_lines, get_content(id_lines, args.batch))
    logger.debug("Begginning to iterate.")
    for i, text_content in enumerate(input_iter):
        batch_list.append(text_content)
        if (i+1) % args.batch is 0:
            outputs.update(read(batch_list, *read_args, **read_kwargs))
            batch_list = []
    logger.debug("Finished iteration.")
    if len(batch_list) > 0:
        logger.debug("Reading remaining files.")
        outputs.update(read(batch_list, *read_args, **read_kwargs))
    with open(path.basename(base_dir) + '_outputs.pkl', 'wb') as f:
        pickle.dump(outputs, f)
