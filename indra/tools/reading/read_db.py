"""This module provides essential tools to run reading using indra's own
database. This may also be run as a script; for details run: 
`python read_pmids_db --help`
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import sys
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
from sqlalchemy.sql.expression import or_, and_, intersect


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
              ' -i/--id_file. For safest use, files should be given by '
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
        help=('Make the reader read all the content. Otherwise, the system '
              'will simply process existing readings where possible.'),
        action='store_true'
        )
    parser.add_argument(
        '--force_fulltext',
        help='Make the reader only read full text from the database.',
        action='store_true'
        )
    parser.add_argument(
        '--no_read',
        help='Only create statements using existing reading results.',
        action='store_true'
        )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

from indra.util import unzip_string, zip_string
from indra.db import get_primary_db, formats, texttypes
from indra.tools.reading.read_pmids import get_mem_total
from indra.util.get_version import get_version as get_indra_version
from indra.sources import reach, sparser


class ReadingError(Exception):
    pass


class ReachError(ReadingError):
    pass


class SparserError(ReadingError):
    pass


# =============================================================================
# Useful functions
# =============================================================================


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


def _enrich_reading_data(reading_data_iter, db=None):
    """Get db ids for all ReadingData objects that correspond to a db ref.

    Note that the objects are modified IN PLACE, so nothing is returned, and if
    a copy of the objects is passed as an argument, this function will have no
    effect. This does nothing if the readings are not in the database.
    """
    logger.debug("Enriching the reading data with database refs.")
    if db is None:
        db = get_primary_db()
    possible_matches = db.select_all(
        'readings',
        db.Readings.text_content_id.in_([rd.tcid for rd in reading_data_iter
                                         if rd.reading_id is None])
        )
    for rdata in reading_data_iter:
        for reading in possible_matches:
            if rdata.matches(reading):
                rdata.reading_id = reading.id
                break
    return


# =============================================================================
# Content Retrieval
# =============================================================================


def get_clauses(id_str_list, db):
    """Get a list of clauses to be passed to a db query."""
    id_types = db.TextRef.__table__.columns.keys()
    id_dict = {id_type: [] for id_type in id_types}
    for id_entry in id_str_list:
        id_type, id_val = _convert_id_entry(id_entry, id_types)
        id_dict[id_type].append(id_val)
    id_condition_list = [getattr(db.TextRef, id_type).in_(id_list)
                         for id_type, id_list in id_dict.items()
                         if len(id_list)]
    return [or_(*id_condition_list)]


def get_text_content_summary_string(q, db):
    """Create a table with some summary data for a query."""
    N_tot = q.count()
    if N_tot > 0:
        log_n = floor(log10(N_tot))
    else:
        log_n = 1
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


def get_content_query(id_str_list, readers, db=None, force_fulltext=False,
                      force_read=False):
    """Load all the content that will be read."""
    if db is None:
        db = get_primary_db()
    logger.debug("Got db handle.")

    tc_tr_binding = db.TextContent.text_ref_id == db.TextRef.id
    rd_tc_binding = db.Readings.text_content_id == db.TextContent.id
    general_clauses = [tc_tr_binding]
    if force_fulltext:
        general_clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)
    id_clauses = get_clauses(id_str_list, db)
    logger.debug("Generated %d id clauses." % len(id_clauses))

    if id_clauses:
        # Get the text content query object
        tc_query = db.filter_query(
            db.TextContent,
            *(general_clauses + id_clauses)
            )
        try:
            logger.debug("Going to try to make a nice summary...")
            logger.info(get_text_content_summary_string(tc_query, db))
        except Exception:
            logger.info("Could not print summary of results.")

        if not force_read:
            logger.debug("Getting content to be read.")
            # Each sub query is a set of content that has been read by one of
            # the readers.
            tc_q_subs = [tc_query.filter(rd_tc_binding, r.matches_clause(db))
                         for r in readers]
            tc_tbr_query = tc_query.except_(intersect(*tc_q_subs))
        else:
            logger.debug('All content will be read (force_read).')
            tc_tbr_query = tc_query
    else:
        logger.info("Did not retreive content or readings from database.")
        tc_tbr_query = None

    return tc_tbr_query.distinct()


def get_readings_query(id_str_list, readers, db=None, force_fulltext=False):
    """Get all the readings available for the id's in id_str_list."""
    if db is None:
        db = get_primary_db()
    general_clauses = [
        # Bind conditions on readings to conditions on content.
        db.Readings.text_content_id == db.TextContent.id,

        # Bind text content to text refs
        db.TextContent.text_ref_id == db.TextRef.id,

        # Check if at least one of the readers has read the content
        or_(*[reader.matches_clause(db) for reader in readers])
        ]
    if force_fulltext:
        general_clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)
    id_clauses = get_clauses(id_str_list, db)
    if id_clauses:
        readings_query = db.filter_query(
            db.Readings,

            # Bind conditions on readings to conditions on content.
            db.Readings.text_content_id == db.TextContent.id,

            # Bind text content to text refs
            db.TextContent.text_ref_id == db.TextRef.id,

            # Check if at least one of the readers has read the content
            or_(*[reader.matches_clause(db) for reader in readers]),

            # Conditions generated from the list of ids. These include a
            # text-ref text-content binding to connect with id data.
            *get_clauses(id_str_list, db)
            )
    else:
        readings_query = None
    return readings_query.distinct()


# =============================================================================
# Reader Classes
# =============================================================================


class Reader(object):
    """This abstract object defines and some general methods for readers."""
    name = NotImplemented

    def __init__(self, base_dir=None, n_proc=1):
        if base_dir is None:
            base_dir = _get_dir('run_' + self.name)
        self.n_proc = n_proc
        self.base_dir = _get_dir(base_dir)
        tmp_dir = tempfile.mkdtemp(
            prefix='reach_job_%s' % _time_stamp(),
            dir=base_dir
            )
        self.tmp_dir = tmp_dir
        self.input_dir = _get_dir(tmp_dir, 'input')
        return

    def read(self, read_list, verbose=False, force_read=True):
        "Read a list of items and return a dict of output files."
        raise NotImplementedError()

    def matches_clause(self, db):
        "Make the clauses to get content that match Reader version and name."
        return and_(db.Readings.reader == self.name,
                    db.Readings.reader_version == self.version)


def get_reader_children():
    """Get all children of the Reader objcet."""
    try:
        children = Reader.__subclasses_()
    except AttributeError:
        module = sys.modules[__name__]
        children = [cls for cls_name, cls in module.__dict__.items()
                    if isinstance(cls, type) and issubclass(cls, Reader)
                    and cls_name != 'Reader']
    return children


class ReachReader(Reader):
    """This object encodes an interface to the reach reading script."""
    REACH_MEM = 5  # GB
    MEM_BUFFER = 2  # GB
    name = 'REACH'

    def __init__(self, *args, **kwargs):
        self.exec_path, self.version = self._check_reach_env()
        super(ReachReader, self).__init__(*args, **kwargs)
        conf_fmt_fname = pjoin(path.dirname(__file__),
                               'reach_conf_fmt.txt')
        self.conf_file_path = pjoin(self.tmp_dir, 'indra.conf')
        with open(conf_fmt_fname, 'r') as fmt_file:
            fmt = fmt_file.read()
            loglevel = 'INFO'  # 'DEBUG' if logger.level == logging.DEBUG else 'INFO'
            with open(self.conf_file_path, 'w') as f:
                f.write(
                    fmt.format(tmp_dir=self.tmp_dir, num_cores=self.n_proc,
                               loglevel=loglevel)
                    )
        self.output_dir = _get_dir(self.tmp_dir, 'output')
        return

    @classmethod
    def _join_json_files(cls, prefix):
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
            # Remove .uaz.<subfile type>.json
            prefix = '.'.join(path.basename(json_file).split('.')[:-3])
            json_prefixes.add(pjoin(self.output_dir, prefix))

        # Join each set of json files and store the json dict.
        reading_data_list = []
        for prefix in json_prefixes:
            base_prefix = path.basename(prefix)
            if base_prefix.isdecimal():
                base_prefix = int(base_prefix)
            reading_data_list.append(ReadingData(
                base_prefix,
                self.name,
                self.version,
                formats.JSON,
                self._join_json_files(prefix)
                ))
            logger.debug('Joined files for prefix %s.' % base_prefix)
        return reading_data_list

    def read(self, read_list, verbose=False):
        """Read the content, returning a dict of ReadingData objects."""
        init_msg = 'Running %s with:\n' % self.name
        init_msg += '\n'.join([
            'n_proc=%s' % self.n_proc
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


class SparserReader(Reader):
    """This object provides methods to interface with the commandline tool."""

    name = 'SPARSER'

    def __init__(self, *args, **kwargs):
        self.version = sparser.get_version()
        super(SparserReader, self).__init__(*args, **kwargs)
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
                # This implies that it is a file path
                fpath = item.strip()
                if fpath.endswith('.nxml'):
                    # If it is already an nxml, we just need to adjust the
                    # name a bit, if anything.
                    if fpath.startswith('PMC'):
                        file_list.append(fpath)
                    else:
                        new_fpath = pjoin(self.tmp_dir, path.basename(fpath))
                        shutil.copy(fpath, new_fpath)
                else:
                    # Otherwise we need to frame the content in xml and put it
                    # in a new file with the appropriat name.
                    old_name = path.basename(fpath)
                    new_fname = '.'.join(old_name.split('.')[:-1] + ['nxml'])
                    new_fpath = pjoin(self.tmp_dir, new_fname)
                    with open(fpath, 'r') as f_old:
                        content = f_old.read()
                    nxml_str = sparser.make_sparser_nxml_from_text(content)
                    with open(new_fpath, 'w') as f_new:
                        f_new.write(nxml_str)
                    file_list.append(new_fpath)
            elif all([hasattr(item, a) for a in ['format', 'content', 'id']]):
                # This implies that it is a text content object, or something
                # with a matching API.
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
        reading_data_list = []
        patt = re.compile(r'(.*?)-semantics.*?')
        for outpath in output_files:
            re_out = patt.match(path.basename(outpath))
            if re_out is None:
                raise SparserError("Could not get prefix from output path %s."
                                   % outpath)
            prefix = re_out.groups()[0]
            if prefix.startswith('PMC'):
                prefix = prefix[3:]
            if prefix.isdecimal():
                # In this case we assume the prefix is a tcid.
                prefix = int(prefix)

            with open(outpath, 'r') as f:
                content = json.load(f)

            reading_data_list.append(ReadingData(
                prefix,
                self.name,
                self.version,
                formats.JSON,
                content
                ))
        return reading_data_list

    def read(self, read_list, verbose=False, log=False):
        "Perform the actual reading."
        ret = None
        file_list = self.prep_input(read_list)
        if len(file_list) > 0:
            logger.info("Beginning to run sparser.")
            output_file_list = []
            if log:
                log_name = 'sparser_run.log'
                outbuf = open(log_name, 'w')
            else:
                outbuf = None
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


# =============================================================================
# Core Reading Functions
# =============================================================================


def read_content(read_list, readers, *args, **kwargs):
    """Perform the reading, returning dicts of jsons."""
    output_list = []
    for reader in readers:
        res_list = reader.read(read_list, *args, **kwargs)
        if res_list is None:
            logger.info("Nothing read by %s." % reader.name)
        else:
            logger.info("Successfully read %d content entries with %s."
                        % (len(res_list), reader.name))
            output_list += res_list
    logger.info("Read %s text content entries in all." % len(output_list))
    return output_list


def read_db(id_str_list, readers, db=None, **kwargs):
    """Read contents retrieved from the database.

    The content will be retrieved in batchs, given by the `batch` argument.
    This prevents the system RAM from being overloaded.

    Parameters
    ----------
    id_str_list : list of stings
        A list of id strings, as would retrieved from a file, each of the form
        '<id type>:<id value>', for example 'pmid:12345'
    readers : list of reader objects
        A list of the readers that will be use, for example ['reach'] if you
        wanted to use the reach reader.
    batch : int
        The number of content entries read for each batch. Default 1000.
    force_fulltext : bool
        If True, only get fulltext content from the database. Default False.

    Other keyword arguments are passed to the `read_content` function.

    Returns
    -------
    outputs : list of ReadingData instances
        The results of the readings with relevant metadata.
    """
    if db is None:
        db = get_primary_db()

    # Retrieve the kwargs needed in this function.
    batch_size = kwargs.pop('batch', 1000)
    force_fulltext = kwargs.pop('force_fulltext', False)
    force_read = kwargs.pop('force_read', False)

    # Get the iterator.
    logger.debug("Getting iterator.")
    tc_read_q = get_content_query(
        id_str_list,
        readers,
        db=db,
        force_fulltext=force_fulltext,
        force_read=force_read
        )
    logger.debug("Begginning to iterate.")
    batch_list_dict = {r.name: [] for r in readers}
    new_outputs = []
    if tc_read_q is not None:
        for text_content in tc_read_q.yield_per(batch_size):
            # The get_content function returns an iterator which yields
            # results in batches, so as not to overwhelm RAM. We need to read
            # in batches for much the same reaason.
            for r in readers:
                # Try to get a previous reading from this reader.
                reading = db.select_one(
                    db.Readings,
                    db.Readings.text_content_id == text_content.id,
                    r.matches_clause(db)
                    )
                if reading is not None and not force_read:
                    continue

                batch_list_dict[r.name].append(text_content)
                if (len(batch_list_dict[r.name])+1) % batch_size is 0:
                    # TODO: this is a bit cludgy...maybe do this better?
                    # Perhaps refactor read_content.
                    logger.debug("Reading batch of files for %s." % r.name)
                    new_outputs += read_content(batch_list_dict[r.name], [r],
                                                **kwargs)
                    batch_list_dict[r.name] = []
        logger.debug("Finished iteration.")
        # Pick up any stragglers.
        for r in readers:
            if len(batch_list_dict[r.name]) > 0:
                logger.debug("Reading remaining files for %s." % r.name)
                new_outputs += read_content(batch_list_dict[r.name], [r],
                                            **kwargs)
    return new_outputs


def get_readings(id_lines, readers, force_fulltext=False, batch_size=1000):
    """Get readings from the database."""
    # Get any previous readings. Note that we do this BEFORE posting the new
    # readings. Otherwise we would have duplicates.
    previous_readings_query = get_readings_query(
        id_lines,
        readers,
        force_fulltext=force_fulltext
        )
    if previous_readings_query is not None:
        prev_readings = [
            ReadingData(
                r.text_content_id,
                r.reader,
                r.reader_version,
                r.format,
                zlib.decompress(r.bytes, 16+zlib.MAX_WBITS).decode('utf8'),
                r.id
                )
            for r in previous_readings_query.yield_per(batch_size)
            ]
    else:
        prev_readings = []
    return prev_readings


def read_files(files, readers, **kwargs):
    """Read the files in `files` with the reader objects in `readers`."""
    return read_content(files, readers, **kwargs)


# =============================================================================
# Reading Processing
# =============================================================================


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

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('text_content_id', 'reader', 'reader_version', 'format',
                'bytes')

    def get_statements(self):
        """General method to create statements."""
        if self.reader == ReachReader.name:
            if self.format == formats.JSON:
                # Process the reach json into statements.
                json_str = json.dumps(self.content)
                stmts = reach.process_json_str(json_str).statements
            elif self.romat == formats.XML:
                # Process the reach xml into statements (untested)
                stmts = reach.process_nxml_str(self.content.to_string())
            else:
                logger.error("Unknown format for creating reach statments: %s."
                             % self.format)
                stmts = []
        elif self.reader == SparserReader.name:
            if self.format == formats.JSON:
                # Process the sparser content into statements
                stmts = sparser.process_json_dict(self.content).statements
            else:
                logger.error("Sparser should only ever be JSON, not %s."
                             % self.format)
                stmts = []
        return stmts

    def zip_content(self):
        """Compress the content, returning bytes."""
        if self.format == formats.JSON:
            ret = zip_string(json.dumps(self.content))
        elif self.format == formats.TEXT:
            ret = zip_string(self.content)
        else:
            raise Exception('Do not know how to zip format %s.' % self.format)
        return ret

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


def post_reading_output(output_list, db=None):
    """Put the reading output on the database."""
    if db is None:
        db = get_primary_db()

    # Create the list of records to be copied, ensuring no uniqueness conflicts
    r_list = db.select_all(
        db.Readings,
        db.Readings.text_content_id.in_([rd.tcid for rd in output_list])
        )
    exisiting_tcid_set = set([r.text_content_id for r in r_list])
    upload_list = []
    for reading_data in output_list:
        # First check if this tcid is even in the set of existing tcids in the
        # readings table.
        if reading_data.tcid in exisiting_tcid_set:
            r_tcid_list = [r for r in r_list
                           if r.text_content_id == reading_data.tcid]
            # Now check for any exact matches:
            if any([reading_data.matches(r) for r in r_tcid_list]):
                continue

        # If there were no conflicts, we can add this to the copy list.
        upload_list.append(reading_data.make_tuple())

    # Copy into the database.
    logger.info("Adding %d/%d reading entries to the database." %
                (len(upload_list), len(output_list)))
    db.copy('readings', upload_list, ReadingData.get_cols())
    return


def produce_readings(input_list, reader_list, input_mode, verbose=False,
                     force_read=False, force_fulltext=False, batch_size=1000,
                     no_read=False, no_upload=False, pickle_result=False):
    """Produce the reading output for the given ids."""
    if no_read:
        return []

    if input_mode is 'ids':
        outputs = read_db(input_list, reader_list, verbose=verbose,
                          force_read=force_read,
                          force_fulltext=force_fulltext, batch=batch_size)
        prev_readings = get_readings(input_list, reader_list, force_fulltext,
                                     batch_size)
    elif input_mode is 'files':
        outputs = read_files(input_list, readers, verbose=verbose)
    else:
        raise ReadingError("Unknown input_mode: %s." % input_mode)

    if pickle_result:
        reading_out_path = pjoin(os.getcwd(), 'reading_outputs.pkl')
        with open(reading_out_path, 'wb') as f:
            pickle.dump([output.make_tuple() for output in outputs], f)
        print("Reading outputs stored in %s." % reading_out_path)

    if not no_upload and input_mode == 'ids':
        post_reading_output(outputs)

    if input_mode is 'ids':
        outputs += prev_readings

    return outputs


# =============================================================================
# Statement Processing
# =============================================================================


class StatementData(object):
    """Contains metadata for statements, as well as the statement itself."""

    def __init__(self, statement, reading_id):
        self.reading_id = reading_id
        self.statement = statement
        self.indra_version = get_indra_version()
        return

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('reader_ref', 'uuid', 'type', 'json', 'indra_version')

    def make_tuple(self):
        """Make a tuple for copying into the database."""
        return (self.reading_id, self.statement.uuid,
                self.statement.__class__.__name__,
                json.dumps(self.statement.to_json()),
                self.indra_version)


def produce_statements(output_list, enrich=True, no_upload=False,
                       pickle_result=False):
    """Convert the reader output into a list of StatementData instances."""
    if enrich:
        _enrich_reading_data(output_list)

    stmt_data_list = [StatementData(stmt, output.reading_id)
                      for output in output_list
                      for stmt in output.get_statements()]
    logger.info("Found %d statements from %d readings." %
                (len(stmt_data_list), len(output_list)))

    if not no_upload:
        upload_statements(stmt_data_list)
    if pickle_result:
        stmts_path = pjoin(os.getcwd(), 'statements.pkl')
        with open(stmts_path, 'wb') as f:
            pickle.dump([sd.statement for sd in stmt_data_list], f)
        print("Statements pickled in %s." % stmts_path)

    return stmt_data_list


def upload_statements(stmt_data_list, db=None):
    """Upload the statements to the database."""
    if db is None:
        db = get_primary_db()

    logger.info("Uploading %d statements to the database." %
                len(stmt_data_list))
    db.copy('statements', [s.make_tuple() for s in stmt_data_list],
            StatementData.get_cols())

    logger.info("Uploading agents to the database.")
    reading_id_set = set([sd.reading_id for sd in stmt_data_list])
    db.insert_agents([sd.statement for sd in stmt_data_list],
                     db.Statements.reader_ref.in_(reading_id_set))
    return


# =============================================================================
# Main for script use
# =============================================================================


if __name__ == "__main__":
    # Process the arguments. =================================================
    if args.id_file and args.file_file:
        raise ReadingError("Cannot process both files and ids.")

    # Get the ids or files.
    if args.id_file is not None:
        with open(args.id_file, 'r') as f:
            input_lines = f.readlines()
        logger.info("Found %d ids." % len(input_lines))
        mode = 'ids'
    elif args.file_file is not None:
        with open(args.file_file, 'r') as f:
            input_lines = f.readlines()
        logger.info("Found %d files." % len(input_lines))
        for ftype in ['nxml', 'txt']:
            logger.debug('%d are %s' % (
                len([f for f in input_lines if f.endswith(ftype)]), ftype
                ))
        mode = 'files'
    else:
        raise ReadingError('No inputs provided.')

    # Select only a sample of the lines, if sample is chosen.
    if args.sample is not None:
        input_lines = random.sample(input_lines, args.sample)

    # If a range is specified, only use that range.
    if args.in_range is not None:
        start_idx, end_idx = [int(n) for n in args.in_range.split(':')]
        input_lines = input_lines[start_idx:end_idx]

    # Create a single base directory
    base_dir = _get_dir('run_%s' % ('_and_'.join(args.readers)))

    # Get the readers objects.
    readers = [reader_class(base_dir=base_dir, n_proc=args.num_procs)
               for reader_class in get_reader_children()
               if reader_class.name.lower() in args.readers]

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    # Read everything ========================================================
    outputs = produce_readings(input_lines, readers, mode, verbose=verbose,
                               force_read=args.force_read,
                               force_fulltext=args.force_fulltext,
                               batch_size=args.batch, no_read=args.no_read,
                               no_upload=args.no_reading_upload,
                               pickle_result=args.pickle)

    # Convert the outputs to statements ======================================
    produce_statements(outputs, enrich=(mode == 'ids'),
                       no_upload=(args.no_statement_upload or mode != 'ids'),
                       pickle_result=args.pickle)
