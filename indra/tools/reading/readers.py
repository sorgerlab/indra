"""Objects for interacting with bulk nlp reading tools."""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import sys
import shutil
import re
import tempfile
import glob
import json
import logging
import subprocess
import zlib
from os import path, mkdir, environ, listdir, remove
from io import BytesIO
from datetime import datetime
from multiprocessing import Pool
from platform import system

from indra.db import formats
from indra.db import sql_expressions as sql
from indra.util import zip_string
from indra.sources import sparser, reach


logger = logging.getLogger('readers')


def _get_dir(*args):
    dirname = path.join(*args)
    if path.isabs(dirname):
        dirpath = dirname
    elif path.exists(dirname):
        dirpath = path.abspath(dirname)
    else:
        dirpath = path.join(path.dirname(path.abspath(__file__)), dirname)
    if not path.exists(dirpath):
        mkdir(dirpath)
    return dirpath


def _time_stamp():
    return datetime.now().strftime("%Y%m%d%H%M%S")


def _get_mem_total():
    if system() == 'Linux':
        with open('/proc/meminfo', 'r') as f:
            lines = f.readlines()
        tot_entry = [line for line in lines if line.startswith('MemTotal')][0]
        ret = int(tot_entry.split(':')[1].replace('kB', '').strip())/10**6
    else:
        ret = None
    return ret


class ReadingError(Exception):
    pass


class ReachError(ReadingError):
    pass


class SparserError(ReadingError):
    pass


class Reader(object):
    """This abstract object defines and some general methods for readers."""
    name = NotImplemented

    def __init__(self, base_dir=None, n_proc=1):
        if base_dir is None:
            base_dir = 'run_' + self.name.lower()
        self.n_proc = n_proc
        self.base_dir = _get_dir(base_dir)
        tmp_dir = tempfile.mkdtemp(
            prefix='%s_job_%s' % (self.name.lower(), _time_stamp()),
            dir=self.base_dir
            )
        self.tmp_dir = tmp_dir
        self.input_dir = _get_dir(tmp_dir, 'input')
        return

    def read(self, read_list, verbose=False, log=False):
        "Read a list of items and return a dict of output files."
        raise NotImplementedError()

    def matches_clause(self, db):
        "Make the clauses to get content that match Reader version and name."
        return sql.and_(db.Readings.reader == self.name,
                        db.Readings.reader_version == self.version[:20])


class ReachReader(Reader):
    """This object encodes an interface to the reach reading script."""
    REACH_MEM = 5  # GB
    MEM_BUFFER = 2  # GB
    name = 'REACH'

    def __init__(self, *args, **kwargs):
        self.exec_path, self.version = self._check_reach_env()
        super(ReachReader, self).__init__(*args, **kwargs)
        conf_fmt_fname = path.join(path.dirname(__file__),
                                   'util/reach_conf_fmt.txt')
        self.conf_file_path = path.join(self.tmp_dir, 'indra.conf')
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
    def _join_json_files(cls, prefix, clear=False):
        """Join different REACH output JSON files into a single JSON object.

        The output of REACH is broken into three files that need to be joined
        before processing. Specifically, there will be three files of the form:
        `<prefix>.uaz.<subcategory>.json`.

        Parameters
        ----------
        prefix : str
            The absolute path up to the extensions that reach will add.
        clear : bool
            Default False - if True, delete the files as soon as they are
            loaded.

        Returns
        -------
        json_obj : dict
            The result of joining the files, keyed by the three subcategories.
        """
        filetype_list = ['entities', 'events', 'sentences']
        json_dict = {}
        try:
            for filetype in filetype_list:
                fname = prefix + '.uaz.' + filetype + '.json'
                with open(fname, 'rt') as f:
                    json_dict[filetype] = json.load(f)
                if clear:
                    remove(fname)
                    logger.debug("Removed %s." % fname)
        except IOError as e:
            logger.error(
                'Failed to open JSON files for %s; REACH error?' % prefix
                )
            logger.exception(e)
            return None
        return json_dict

    def _check_reach_env(self):
        """Check that the environment supports runnig reach."""
        # Get the path to the REACH JAR
        path_to_reach = environ.get('REACHPATH', None)
        if path_to_reach is None or not path.exists(path_to_reach):
            raise ReachError(
                'Reach path unset or invalid. Check REACHPATH environment var.'
                )

        logger.debug('Using REACH jar at: %s' % path_to_reach)

        # Get the reach version.
        reach_version = environ.get('REACH_VERSION', None)
        if reach_version is None:
            logger.debug('REACH version not set in REACH_VERSION')
            m = re.match('reach-(.*?)\.jar', path.basename(path_to_reach))
            reach_version = re.sub('-SNAP.*?$', '', m.groups()[0])

        logger.debug('Using REACH version: %s' % reach_version)
        return path_to_reach, reach_version

    def write_content(self, text_content):
        def write_content_file(ext):
            fname = '%s.%s' % (text_content.id, ext)
            with open(path.join(self.input_dir, fname), 'wb') as f:
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
                    path.join(self.input_dir, path.basename(fname))
                    )
            else:
                self.write_content(text_content)
        return

    def get_output(self):
        """Get the output of a reading job as a list of filenames."""
        logger.info("Getting outputs.")
        # Get the set of prefixes (each will correspond to three json files.)
        json_files = glob.glob(path.join(self.output_dir, '*.json'))
        json_prefixes = set()
        for json_file in json_files:
            # Remove .uaz.<subfile type>.json
            prefix = '.'.join(path.basename(json_file).split('.')[:-3])
            json_prefixes.add(path.join(self.output_dir, prefix))

        # Join each set of json files and store the json dict.
        reading_data_list = []
        for prefix in json_prefixes:
            base_prefix = path.basename(prefix)
            if base_prefix.isdecimal():
                base_prefix = int(base_prefix)
            try:
                content = self._join_json_files(prefix, clear=True)
            except Exception as e:
                logger.exception(e)
                logger.error("Coule not load result for prefix %s." % prefix)
                continue
            reading_data_list.append(ReadingData(
                base_prefix,
                self.name,
                self.version,
                formats.JSON,
                content
                ))
            logger.debug('Joined files for prefix %s.' % base_prefix)
        return reading_data_list

    def clear_input(self):
        """Remove all the input files (at the end of a reading)."""
        for item in listdir(self.input_dir):
            item_path = path.join(self.input_dir, item)
            if path.isfile(item_path):
                remove(item_path)
                logger.debug('Removed input %s.' % item_path)
        return

    def read(self, read_list, verbose=False, log=False):
        """Read the content, returning a list of ReadingData objects."""
        ret = []
        mem_tot = _get_mem_total()
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
            log_file_str = ''
            for line in iter(p.stdout.readline, b''):
                log_line = 'REACH: ' + line.strip().decode('utf8')
                if verbose:
                    logger.info(log_line)
                if log:
                    log_file_str += log_line + '\n'
            if log:
                with open('reach_run.log', 'ab') as f:
                    f.write(log_file_str.encode('utf8'))
            p_out, p_err = p.communicate()
            if p.returncode:
                logger.error('Problem running REACH:')
                logger.error('Stdout: %s' % p_out.decode('utf-8'))
                logger.error('Stderr: %s' % p_err.decode('utf-8'))
                raise ReachError("Problem running REACH")
            logger.info("Reach finished.")
            ret = self.get_output()
            self.clear_input()
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
            fpath = path.join(self.input_dir, 'PMC%d.nxml' % tcid)
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
                        new_fpath = path.join(self.tmp_dir,
                                              'PMC' + path.basename(fpath))
                        shutil.copy(fpath, new_fpath)
                        file_list.append(new_fpath)
                else:
                    # Otherwise we need to frame the content in xml and put it
                    # in a new file with the appropriat name.
                    old_name = path.basename(fpath)
                    new_fname = '.'.join(old_name.split('.')[:-1] + ['nxml'])
                    new_fpath = path.join(self.tmp_dir, new_fname)
                    with open(fpath, 'r') as f_old:
                        content = f_old.read()
                    nxml_str = sparser.make_nxml_from_text(content)
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
                    nxml_str = sparser.make_nxml_from_text(
                        txt_bts.decode('utf8')
                        )
                    add_nxml_file(item.id, nxml_str.encode('utf8'))
                else:
                    raise SparserError("Unrecognized format %s." % item.format)
            else:
                raise SparserError("Unknown type of item for reading %s." %
                                   type(item))
        return file_list

    def get_output(self, output_files, clear=True):
        "Get the output files as an id indexed dict."
        reading_data_list = []
        patt = re.compile(r'(.*?)-semantics.*?')
        for outpath in output_files:
            if outpath is None:
                logger.warning("Found outpath with value None. Skipping.")
                continue

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

            try:
                with open(outpath, 'rt') as f:
                    content = json.load(f)
            except Exception as e:
                logger.exception(e)
                logger.error("Could not load reading content from %s."
                             % outpath)
                continue

            reading_data_list.append(ReadingData(
                prefix,
                self.name,
                self.version,
                formats.JSON,
                content
                ))
            if clear:
                try:
                    remove(outpath)
                    input_path = outpath.replace('-semantics.json', '.nxml')
                    remove(input_path)
                except Exception as e:
                    logger.exception(e)
                    logger.error("Could not remove sparser files %s and %s."
                                 % (outpath, input_path))
        return reading_data_list

    def read_one(self, fpath, outbuf=None, verbose=False):
        if outbuf is None:
            outbuf = BytesIO()
        outbuf.write(b'\nReading %s.\n' % fpath.encode('utf8'))
        outbuf.flush()
        if verbose:
            logger.info('Reading %s.' % fpath)
        outpath = None
        try:
            outpath = sparser.run_sparser(fpath, 'json', outbuf)
        except Exception as e:
            if verbose:
                logger.error('Failed to run sparser on %s.' %
                             fpath)
                logger.exception(e)
            outbuf.write(b'Reading failed.\n')
        return outpath, outbuf

    def read_some(self, fpath_list, outbuf=None, verbose=False):
        "Perform a few readings."
        outpath_list = []
        for fpath in fpath_list:
            output, outbuf = self.read_one(fpath, outbuf, verbose)
            if output is not None:
                outpath_list.append(output)
        return outpath_list, outbuf

    def read(self, read_list, verbose=False, log=False, n_per_proc=None):
        "Perform the actual reading."
        ret = []
        file_list = self.prep_input(read_list)
        if len(file_list) > 0:
            logger.info("Beginning to run sparser.")
            output_file_list = []
            if log:
                log_name = 'sparser_run_%s.log' % _time_stamp()
                outbuf = open(log_name, 'wb')
            else:
                outbuf = None
            try:
                if self.n_proc == 1:
                    for fpath in file_list:
                        outpath, _ = self.read_one(fpath, outbuf, verbose)
                        if outpath is not None:
                            output_file_list.append(outpath)
                else:
                    if n_per_proc is None:
                        L = len(file_list)
                        n_per_proc = max(1, min(1000, L//self.n_proc//2))

                    try:
                        pool = Pool(self.n_proc)
                        if n_per_proc is not 1:
                            batches = [file_list[n*n_per_proc:(n+1)*n_per_proc]
                                       for n in range(L//n_per_proc + 1)]
                            out_lists_and_buffs = pool.map(self.read_some,
                                                           batches)
                        else:
                            out_files_and_buffs = pool.map(self.read_one,
                                                           file_list)
                            out_lists_and_buffs = [([out_files], buffs)
                                                   for out_files, buffs
                                                   in out_files_and_buffs]
                    finally:
                        pool.close()
                        pool.join()
                    for i, (out_list, buff) in enumerate(out_lists_and_buffs):
                        if out_list is not None:
                            output_file_list += out_list
                        if log:
                            outbuf.write(b'Log for producing output %d/%d.\n'
                                         % (i, len(out_lists_and_buffs)))
                            buff.seek(0)
                            outbuf.write(buff.read() + b'\n')
                            outbuf.flush()
            finally:
                if log:
                    outbuf.close()
                    if verbose:
                        logger.info("Sparser logs may be found at %s." %
                                    log_name)
            ret = self.get_output(output_file_list)
        return ret


def get_readers():
    """Get all children of the Reader objcet."""
    try:
        children = Reader.__subclasses__()
    except AttributeError:
        module = sys.modules[__name__]
        children = [cls for cls_name, cls in module.__dict__.items()
                    if isinstance(cls, type) and issubclass(cls, Reader)
                    and cls_name != 'Reader']
    return children


class ReadingData(object):
    """Object to contain the data produced by a reading.

    This is primarily designed for use with the database.

    Init Parameters
    ---------------
    tcid : int or str
        An identifier of the text content that produced the reading. Must
        be an int for use with the database.
    reader : str
        The name of the reader, consistent with it's `name` attribute, for
        example: 'REACH'
    reader_version : str
        A string identifying the version of the underlying nlp reader.
    content_format : str
        The format of the content. Options are in indra.db.formats.
    content : str
        The content of the reading result. A string in the format given by
        `content_format`.
    reading_id : int or None
        Optional. The id corresponding to the Readings entry in the db.
    """

    def __init__(self, tcid, reader, reader_version, content_format, content,
                 reading_id=None):
        self.reading_id = reading_id
        self.tcid = tcid
        self.reader = reader
        self.reader_version = reader_version
        self.format = content_format
        self.content = content
        return

    @classmethod
    def get_cols(self):
        """Get the columns for the tuple returned by `make_tuple`."""
        return ('text_content_id', 'reader', 'reader_version', 'format',
                'bytes')

    def get_statements(self):
        """General method to create statements."""
        logger.debug("Making statements from %s." % self.reading_id)
        if self.reader == ReachReader.name:
            if self.format == formats.JSON:
                # Process the reach json into statements.
                json_str = json.dumps(self.content)
                processor = reach.process_json_str(json_str)
            else:
                raise ReadingError("Incorrect format for Reach output: %s."
                                   % self.format)
        elif self.reader == SparserReader.name:
            if self.format == formats.JSON:
                # Process the sparser content into statements
                processor = sparser.process_json_dict(self.content)
                if processor is not None:
                    processor.set_statements_pmid(None)
            else:
                raise ReadingError("Sparser should only ever be JSON, not %s."
                                   % self.format)
        if processor is None:
            logger.error("Production of statements from %s failed for %s."
                         % (self.reader, self.tcid))
            stmts = []
        else:
            stmts = processor.statements
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
        # Note the temporary fix in clipping the reader version length. This is
        # because the version is for some reason clipped in the database.
        return (r_entry.text_content_id == self.tcid
                and r_entry.reader == self.reader
                and r_entry.reader_version == self.reader_version[:20])
