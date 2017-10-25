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

logger = logging.getLogger('read_db')
if __name__ == '__main__':
    parser = ArgumentParser(
        description='A tool to read and process content from the database.'
        )
    parser.add_argument(
        '-i', '--id_file',
        help=('A file containing a list of ids of the form <id_type>:<id>.'
              'Cannot be used in conjunction with -f/--file_file.')
        )
    parser.add_argument(
        '-f', '--file_file',
        help=('A file containing a list of files to be input into reach. These'
              'should be nxml or txt files. Cannot be used in conjunction with'
              '-i/--id_file. For safesed use, files should be given by'
              'absolute paths.')
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
        '-n', '--num_procs',
        help='Select the number of processes to use.',
        type=int,
        default=1
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
    parser.add_argument(
        '-m', '--messy',
        help='Do not clean up directories created while reading.',
        action='store_true'
        )
    parser.add_argument(
        '-p', '--pickle',
        help='Pickle all results and save in .pkl files.',
        action='store_true'
        )
    parser.add_argument(
        '--no_reading_upload',
        help='Choose not to upload the reading output to the database.',
        action='store_true'
        )
    parser.add_argument(
        '--no_statement_upload',
        help='Choose not to upload the statements to the databse.',
        action='store_true'
        )
    parser.add_argument(
        '--force_read',
        help='Make the reader read all the content. (Not yet implemented.)',
        action='store_true'
        )
    parser.add_argument(
        '--force_fulltext',
        help='Make the reader only read full text from the database.',
        action='store_true'
        )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

from indra.util import unzip_string, zip_string
from indra.db import get_primary_db, formats, texttypes
from indra.tools.reading.read_pmids import get_mem_total
from indra.sources import reach, sparser


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


def join_json_files(prefix):
    """Join different REACH output JSON files into a single JSON object.

    The output of REACH is broken into three files that need to be joined
    before processing. Specifically, there will be three files of the form:
    `<prefix>.uaz.<subcategory>.json`.

    Parameters
    ----------
    prefix : str
        The absolute path up to the extensions that reach will add.

    Returns
    -------
    json_obj : dict
        The result of joining the files, keyed by the three subcategories.
    """
    try:
        with open(prefix + '.uaz.entities.json', 'rt') as f:
            entities = json.load(f)
        with open(prefix + '.uaz.events.json', 'rt') as f:
            events = json.load(f)
        with open(prefix + '.uaz.sentences.json', 'rt') as f:
            sentences = json.load(f)
    except IOError as e:
        logger.error(
            'Failed to open JSON files for %s; REACH error?' % prefix
            )
        logger.exception(e)
        return None
    return {'events': events, 'entities': entities, 'sentences': sentences}


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


def get_content(id_str_list, batch_size=1000, db=None, force_fulltext=False):
    """Load all the content that will be read."""
    if db is None:
        db = get_primary_db()
    logger.debug("Got db handle.")
    clauses = get_clauses(id_str_list, db.TextRef)
    if force_fulltext:
        clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)
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


class ReadingData(object):
    """Object to contain the data produced by a reading."""

    def __init__(self, tcid, reader, reader_version, output_format, content,
                 reading_id=None):
        self.reading_id = reading_id
        self.tcid = tcid
        self.reader = reader
        self.reader_version = reader_version
        self.format = output_format
        self.content = content
        return

    def zip_content(self):
        """Compress the content, returning bytes."""
        if self.format == formats.JSON:
            ret = zip_string(json.dumps(self.content))
        elif self.format == formats.TEXT:
            ret = zip_string(self.content)
        else:
            raise Exception('Do not know how to zip format %s.' % self.format)
        return ret

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('text_content_id', 'reader', 'reader_version', 'format',
                'bytes')

    def make_tuple(self):
        """Make the tuple expected by the database."""
        return (self.tcid, self.reader, self.reader_version, self.format,
                self.zip_content())

    def matches(self, r_entry):
        """Determine if reading data matches the a reading entry from the db.

        Returns True if tcid, reader, reader_version match the corresponding
        elements of a db.Reading instance, else False.
        """
        return (r_entry.text_content_id == self.tcid
                and r_entry.reader == self.reader
                and r_entry.reader_version == self.reader_version)


class ReachError(Exception):
    pass


class ReachReader(object):
    """This object encodes an interface to the reach reading script."""
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
            prefix='reach_job_%s' % _time_stamp(),
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
            tc_id_fname_dict[int(base_prefix)] = ReadingData(
                int(base_prefix),
                self.name,
                self.version,
                formats.JSON,
                join_json_files(prefix)
                )
            logger.debug('Joined files for prefix %s.' % base_prefix)
        return tc_id_fname_dict

    def read(self, read_list, verbose=False, force_read=True):
        """Read the content, returning a dict of ReadingData objects."""
        init_msg = 'Running %s with:\n' % self.name
        init_msg += '\n'.join([
            'n_proc=%s' % self.n_proc,
            'force_read=%s' % force_read
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


class SparserError(Exception):
    pass


class SparserReader(object):
    """This object provides methods to interface with the commandline tool."""

    name = 'SPARSER'

    def __init__(self, base_dir=None, n_proc=1):
        if base_dir is None:
            base_dir = 'run_' + self.name
        self.n_proc = n_proc
        self.base_dir = _get_dir(base_dir)
        self.version = sparser.get_version()
        tmp_dir = tempfile.mkdtemp(
            prefix='sparser_job_%s' % _time_stamp(),
            dir=base_dir
            )
        self.tmp_dir = tmp_dir
        self.input_dir = _get_dir(tmp_dir, 'input')
        return

    def prep_input(self, read_list):
        "Prepare the list of files or text content objects to be read."
        logger.info('Prepping input for sparser.')

        file_list = []

        def add_nxml_file(tcid, nxml_bts):
            fpath = pjoin(self.input_dir, 'PMC%d.nxml' % tcid)
            with open(fpath, 'wb') as f:
                f.write(nxml_bts)
            file_list.append(fpath)

        for item in read_list:
            if isinstance(item, str):
                raise SparserError(
                    "This feature not yet implemented for sparser."
                    )  # TODO: Implement this use-case.
            elif all([hasattr(item, a) for a in ['format', 'content', 'id']]):
                if item.format == formats.XML:
                    add_nxml_file(
                        item.id,
                        zlib.decompress(item.content, 16+zlib.MAX_WBITS)
                        )
                elif item.format == formats.TEXT:
                    txt_bts = zlib.decompress(item.content, 16+zlib.MAX_WBITS)
                    nxml_str = sparser.make_sparser_nxml_from_text(
                        txt_bts.decode('utf8')
                        )
                    add_nxml_file(item.id, nxml_str.encode('utf8'))
                elif item.format == formats.JSON:
                    raise SparserError("I don't know how to handle JSON.")
                else:
                    raise SparserError("Unrecognized format %s." % item.format)
            else:
                raise SparserError("Unknown type of item for reading %s." %
                                   type(item))
        return file_list

    def get_output(self, output_files):
        "Get the output files as an id indexed dict."
        tcid_fpath_dict = {}
        patt = re.compile(r'PMC(\w+)-semantics.*?')
        for outpath in output_files:
            re_out = patt.match(path.basename('outpath'))
            if re_out is None:
                raise SparserError("Could not get id from output path %s." %
                                   outpath)
            tcid = re_out.groups()[0]
            tcid_fpath_dict[tcid] = outpath
        return tcid_fpath_dict

    def read(self, read_list, verbose=False, force_read=True, log=False):
        "Perform the actual reading."
        ret = None
        file_list = self.prep_input(read_list)
        if len(read_list) > 1:
            logger.info("Beginning to run sparser.")
            output_file_list = []
            if log:
                log_name = 'sparser_run.log'
                outbuf = open(log_name, 'w')
            try:
                for fpath in file_list:
                    if log:
                        outbuf.write('\nReading %s.\n' % fpath)
                        outbuf.flush()
                    if verbose:
                        logger.info('Reading %s.' % fpath)
                    try:
                        outpath = sparser.run_sparser(fpath, 'json', outbuf)
                        output_file_list.append(outpath)
                    except Exception as e:
                        if verbose:
                            logger.error('Failed to run sparser on %s.' %
                                         fpath)
                            logger.exception(e)
                        if log:
                            outbuf.write('Reading failed.\n')
            finally:
                if log:
                    outbuf.close()
                    if verbose:
                        logger.info("Sparser logs may be found at %s." %
                                    log_name)
            ret = self.get_output(output_file_list)
        return ret


def read_content(read_list, readers, *args, **kwargs):
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


def read_db(id_str_list, readers, **kwargs):
    """Read contents retrieved from the database.

    The content will be retrieved in batchs, given by the `batch` argument.
    This prevents the system RAM from being overloaded.

    Parameters
    ----------
    id_str_list : list of stings
        A list of id strings, as would retrieved from a file, each of the form
        '<id type>:<id value>', for example 'pmid:12345'
    readers : list of strings
        A list of the readers that will be use, for example ['reach'] if you
        wanted to use the reach reader.
    batch : int
        The number of content entries read for each batch. Default 1000.
    force_fulltext : bool
        If True, only get fulltext content from the database. Default False.

    Other keyword arguments are passed to the `read_content` function.

    Returns
    -------
    outputs : dict of ReadingData instances
        The results of the readings with relevant metadata.
    """
    batch_size = kwargs.pop('batch', 1000)
    force_fulltext = kwargs.pop('force_fulltext', False)
    logger.debug("Getting iterator.")
    input_iter = get_content(id_str_list, batch_size=batch_size,
                             force_fulltext=force_fulltext)
    logger.debug("Begginning to iterate.")
    batch_list = []
    outputs = {}
    for i, text_content in enumerate(input_iter):
        # The get_content function returns an iterator which yields results in
        # batches, so as not to overwhelm RAM. We need to read in batches for
        # much the same reaason.
        batch_list.append(text_content)
        if (i+1) % args.batch is 0:
            outputs.update(read_content(batch_list, readers, **kwargs))
            batch_list = []
    logger.debug("Finished iteration.")
    # Pick up any stragglers.
    if len(batch_list) > 0:
        logger.debug("Reading remaining files.")
        outputs.update(read_content(batch_list, readers, **kwargs))
    return outputs


def read_files(file_str_list, readers, **kwargs):
    """Read the files provided by the list of files."""
    return read_content(file_str_list, readers, **kwargs)


def enrich_reading_data(reading_data_iter, db=None):
    """Get db ids for all ReadingData objects that correspond to a db ref.

    Note that the objects are modified IN PLACE, so nothing is returned, and if
    a copy of the objects is passed as an argument, this function will have no
    effect.
    """
    logging.debug("Enriching the reading data with database refs.")
    if db is None:
        db = get_primary_db()
    possible_matches = db.select_all(
        'readings',
        db.Readings.text_content_id.in_([rd.tcid for rd in reading_data_iter])
        )
    for rdata in reading_data_iter:
        for reading in possible_matches:
            if rdata.matches(reading):
                rdata.reading_id = reading.id
                break
    return


def post_reading_output(output_dict, db=None):
    """Put the reading output on the database."""
    if db is None:
        db = get_primary_db()

    # Create the list of records to be copied, ensuring no uniqueness conflicts
    r_list = db.select_all(
        db.Readings,
        db.Readings.text_content_id.in_(list(output_dict.keys()))
        )
    exisiting_tcid_set = set([r.text_content_id for r in r_list])
    upload_list = []
    for tcid, reading_data in output_dict.items():
        # First check if this tcid is even in the set of existing tcids in the
        # readings table.
        if tcid in exisiting_tcid_set:
            r_tcid_list = [r for r in r_list if r.text_content_id == tcid]
            # Now check for any exact matches:
            if any([reading_data.matches(r) for r in r_tcid_list]):
                continue

        # If there were no conflicts, we can add this to the copy list.
        upload_list.append(reading_data.make_tuple())

    # Copy into the database.
    logging.info("Adding %d/%d reading entries to the database." %
                 (len(upload_list), len(output_dict)))
    db.copy('readings', upload_list, ReadingData.get_cols())
    return


class StatementData(object):
    """Contains metadata for statements, as well as the statement itself."""
    def __init__(self, statement, reading_data):
        self.reading_data = reading_data
        self.statement = statement
        return

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('reader_ref', 'uuid', 'type', 'json')

    def make_tuple(self):
        """Make a tuple for copying into the database."""
        assert self.reading_data.reading_id is not None, \
            "Reading data must be loaded into the database first."
        return (self.reading_data.reading_id, self.statement.uuid,
                self.statement.__class__.__name__,
                json.dumps(self.statement.to_json()))


def make_statements(output_dict):
    """Convert the reader output into a list of StatementData instances."""
    stmts = []
    enrich_reading_data(output_dict.values())
    for output in output_dict.values():
        if output.reader == ReachReader.name:
            # Convert the JSON object into a string first so that a series of
            # string replacements can happen in the REACH processor
            reach_proc = reach.process_json_str(json.dumps(output.content))
            stmts += [StatementData(stmt, output)
                      for stmt in reach_proc.statements]
    logger.info("Found %d statements from %d readings." %
                (len(stmts), len(output_dict)))
    return stmts


def upload_statements(stmts, db=None):
    """Upload the statements to the database."""
    if db is None:
        db = get_primary_db()
    logging.info("Uploading %d statements to the database." % len(stmts))
    db.copy('statements', [s.make_tuple() for s in stmts],
            StatementData.get_cols())
    return


if __name__ == "__main__":
    # Process the arguments. =================================================
    assert not (args.id_file and args.file_file), \
        "Cannot process both files and ids."
    # Get the ids.
    if args.id_file is not None:
        with open(args.id_file, 'r') as f:
            id_lines = f.readlines()
        logger.info("Found %d ids to read." % len(id_lines))
        mode = 'ids'
    elif args.file_file is not None:
        with open(args.file_file, 'r') as f:
            file_lines = f.readlines()
        logger.info("Found %d files to read." % len(file_lines))
        for ftype in ['nxml', 'txt']:
            logger.debug('%d are %s' % (
                len([f for f in file_lines if f.endswith(ftype)]), ftype
                ))
        mode = 'files'
    else:
        raise Exception('No inputs provided.')

    # Select only a sample of the lines, if sample is chosen.
    if args.sample is not None and mode == 'ids':
        id_lines = random.sample(id_lines, args.sample)
        # TODO: Figure out how to handle this in conjunction nwith a file list.

    # If a range is specified, only use that range.
    if args.in_range is not None and mode == 'ids':
        start_idx, end_idx = [int(n) for n in args.in_range.split(':')]
        id_lines = id_lines[start_idx:end_idx]

    # Read everything ========================================================
    if len(id_lines):
        outputs = read_db(id_lines, args.readers,
                          verbose=args.verbose and not args.quiet,
                          force_read=args.force_read,
                          force_fulltext=args.force_fulltext,
                          batch=args.batch, n_proc=args.num_procs)
    else:
        outputs = read_files(file_lines, args.readers,
                             verbose=args.verbose and not args.quiet,
                             n_proc=args.num_procs)

    if args.pickle:
        with open(os.getcwd() + 'reading_outputs.pkl', 'wb') as f:
            pickle.dump(outputs, f)

    if not args.no_reading_upload and mode == 'ids':
        post_reading_output(outputs)

    # Convert the outputs to statements ======================================
    stmts = make_statements(outputs)
    if not args.no_statement_upload and mode == 'ids':
        upload_statements(stmts)
    if args.pickle:
        with open(os.getcwd() + 'statements.pkl', 'wb') as f:
            pickle.dump(stmts, f)
