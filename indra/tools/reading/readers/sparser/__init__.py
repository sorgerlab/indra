import re
import json
import logging

from io import BytesIO
from os import path, remove
from multiprocessing import Pool

from indra.tools.reading.readers.core import Reader, ReadingError
from indra.tools.reading.readers.content import Content
from indra.tools.reading.readers.util import get_time_stamp

from indra.sources import sparser

logger = logging.getLogger(__name__)


class SparserError(ReadingError):
    pass


class SparserReader(Reader):
    """This object provides methods to interface with the commandline tool."""

    name = 'SPARSER'

    def __init__(self, *args, **kwargs):
        self.version = self.get_version()
        super(SparserReader, self).__init__(*args, **kwargs)
        self.file_list = None
        return

    @classmethod
    def get_version(cls):
        return sparser.get_version()

    def prep_input(self, read_list):
        "Prepare the list of files or text content objects to be read."
        logger.info('Prepping input for sparser.')

        self.file_list = []

        for content in read_list:
            quality_issue = self._check_content(content.get_text())
            if quality_issue is not None:
                logger.warning("Skipping %d due to: %s"
                               % (content.get_id(), quality_issue))
                continue

            if content.is_format('nxml'):
                # If it is already an nxml, we just need to adjust the
                # name a bit, if anything.
                if not content.get_filename().startswith('PMC'):
                    content.change_id('PMC' + str(content.get_id()))
                fpath = content.copy_to(self.tmp_dir)
                self.file_list.append(fpath)
            elif content.is_format('txt', 'text'):
                # Otherwise we need to frame the content in xml and put it
                # in a new file with the appropriate name.
                nxml_str = sparser.make_nxml_from_text(content.get_text())
                new_content = Content.from_string('PMC' + str(content.get_id()),
                                                  'nxml', nxml_str)
                fpath = new_content.copy_to(self.tmp_dir)
                self.file_list.append(fpath)
            else:
                raise SparserError("Unrecognized format %s."
                                   % content.format)
        return

    def _map_id(self, content_id):
        if content_id.startswith('PMC'):
            content_id = content_id[3:]
        content_id = super(SparserReader, self)._map_id(content_id)
        return content_id

    def get_output(self, output_files, clear=True):
        "Get the output files as an id indexed dict."
        patt = re.compile(r'(.*?)-semantics.*?')
        for outpath in output_files:
            # Get the reading results if possible.
            if outpath is None:
                logger.warning("Found outpath with value None. Skipping.")
                continue

            try:
                with open(outpath, 'rt') as f:
                    reading = json.load(f)
            except Exception as e:
                logger.exception(e)
                logger.error("Could not load reading content from %s."
                             % outpath)
                reading = None

            # Get the content ID
            re_out = patt.match(path.basename(outpath))
            if re_out is None:
                raise SparserError("Could not get prefix from output path %s."
                                   % outpath)

            content_id = re_out.groups()[0]
            self.add_result(content_id, reading)

            # Clean up the input and output files.
            if clear:
                input_path = outpath.replace('-semantics.json', '.nxml')
                try:
                    remove(outpath)
                    remove(input_path)
                except Exception as e:
                    logger.exception(e)
                    logger.error("Could not remove sparser files %s and %s."
                                 % (outpath, input_path))
        return self.results

    def read_one(self, fpath, outbuf=None, verbose=False):
        fpath = path.abspath(fpath)
        if outbuf is None:
            outbuf = BytesIO()
        outbuf.write(b'\nReading %s.\n' % fpath.encode('utf8'))
        outbuf.flush()
        if verbose:
            logger.info('Reading %s.' % fpath)
        outpath = None
        try:
            outpath = sparser.run_sparser(fpath, 'json', outbuf, timeout=60)
        except Exception as e:
            if verbose:
                logger.error('Failed to run sparser on %s.' %
                             fpath)
                logger.exception(e)
            outbuf.write(b'Reading failed.----------\n')
            outbuf.write(str(e).encode('utf-8') + b'\n')
            outbuf.write(b'-------------------------\n')
        return outpath, outbuf

    def read_some(self, fpath_list, outbuf=None, verbose=False):
        "Perform a few readings."
        outpath_list = []
        for fpath in fpath_list:
            output, outbuf = self.read_one(fpath, outbuf, verbose)
            if output is not None:
                outpath_list.append(output)
        return outpath_list, outbuf

    def _read(self, read_list, verbose=False, log=False, n_per_proc=None):
        "Perform the actual reading."
        ret = []
        self.prep_input(read_list)
        L = len(self.file_list)
        if L == 0:
            return ret

        logger.info("Beginning to run sparser.")
        output_file_list = []
        if log:
            log_name = 'sparser_run_%s.log' % get_time_stamp()
            outbuf = open(log_name, 'wb')
        else:
            outbuf = None
        try:
            if self.n_proc == 1:
                for fpath in self.file_list:
                    outpath, _ = self.read_one(fpath, outbuf, verbose)
                    if outpath is not None:
                        output_file_list.append(outpath)
            else:
                if n_per_proc is None:
                    n_per_proc = max(1, min(1000, L//self.n_proc//2))
                pool = None
                try:
                    pool = Pool(self.n_proc)
                    if n_per_proc is not 1:
                        batches = [self.file_list[n*n_per_proc:(n+1)*n_per_proc]
                                   for n in range(L//n_per_proc + 1)]
                        out_lists_and_buffs = pool.map(self.read_some,
                                                       batches)
                    else:
                        out_files_and_buffs = pool.map(self.read_one,
                                                       self.file_list)
                        out_lists_and_buffs = [([out_files], buffs)
                                               for out_files, buffs
                                               in out_files_and_buffs]
                finally:
                    if pool is not None:
                        pool.close()
                        pool.join()
                for i, (out_list, buff) in enumerate(out_lists_and_buffs):
                    if out_list is not None:
                        output_file_list += out_list
                    if log:
                        outbuf.write(b'Log for producing output %d/%d.\n'
                                     % (i, len(out_lists_and_buffs)))
                        if buff is not None:
                            buff.seek(0)
                            outbuf.write(buff.read() + b'\n')
                        else:
                            outbuf.write(b'ERROR: no buffer was None. '
                                         b'No logs available.\n')
                        outbuf.flush()
        finally:
            if log:
                outbuf.close()
                if verbose:
                    logger.info("Sparser logs may be found at %s." %
                                log_name)
        ret = self.get_output(output_file_list)
        return ret

    @staticmethod
    def get_processor(content):
        processor = sparser.process_json_dict(content)
        if processor is not None:
            processor.set_statements_pmid(None)
        return processor
