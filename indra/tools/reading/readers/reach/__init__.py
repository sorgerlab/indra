import re
import json
import glob
import logging
import subprocess

from os import path, remove, environ, listdir

from indra.config import get_config
from indra.tools.reading.readers.util import get_dir, get_mem_total
from indra.tools.reading.readers.core import Reader, ReadingError

from indra.sources import reach

logger = logging.getLogger(__name__)


class ReachError(ReadingError):
    pass


class ReachReader(Reader):
    """This object encodes an interface to the reach reading script."""
    REACH_MEM = 5  # GB
    MEM_BUFFER = 2  # GB
    name = 'REACH'

    def __init__(self, *args, **kwargs):
        self.exec_path, self.version = self._check_reach_env()
        super(ReachReader, self).__init__(*args, **kwargs)
        conf_fmt_fname = path.join(path.dirname(__file__),
                                   'reach_conf_fmt.txt')
        self.conf_file_path = path.join(self.tmp_dir, 'indra.conf')
        with open(conf_fmt_fname, 'r') as fmt_file:
            fmt = fmt_file.read()
            log_level = 'INFO'
            # log_level = 'DEBUG' if logger.level == logging.DEBUG else 'INFO'
            with open(self.conf_file_path, 'w') as f:
                f.write(
                    fmt.format(tmp_dir=self.tmp_dir, num_cores=self.n_proc,
                               loglevel=log_level)
                )
        self.output_dir = get_dir(self.tmp_dir, 'output')
        self.num_input = 0
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

    @staticmethod
    def _check_reach_env():
        """Check that the environment supports runnig reach."""
        # Get the path to the REACH JAR
        path_to_reach = get_config('REACHPATH')
        if path_to_reach is None:
            path_to_reach = environ.get('REACHPATH', None)
        if path_to_reach is None or not path.exists(path_to_reach):
            raise ReachError(
                'Reach path unset or invalid. Check REACHPATH environment var '
                'and/or config file.'
            )

        logger.debug('Using REACH jar at: %s' % path_to_reach)

        # Get the reach version.
        reach_version = get_config('REACH_VERSION')
        if reach_version is None:
            reach_version = environ.get('REACH_VERSION', None)
        if reach_version is None:
            logger.debug('REACH version not set in REACH_VERSION')
            m = re.match('reach-(.*?)\.jar', path.basename(path_to_reach))
            reach_version = re.sub('-SNAP.*?$', '', m.groups()[0])

        logger.debug('Using REACH version: %s' % reach_version)
        return path_to_reach, reach_version

    @classmethod
    def get_version(cls):
        _, version = cls._check_reach_env()
        return version

    def prep_input(self, read_list):
        """Apply the readers to the content."""
        logger.info("Prepping input.")
        i = 0
        for content in read_list:
            # Check the quality of the text, and skip if there are any issues.
            quality_issue = self._check_content(content.get_text())
            if quality_issue is not None:
                logger.warning("Skipping %d due to: %s"
                               % (content.get_id(), quality_issue))
                continue

            # Look for things that are more like file names, rather than ids.
            cid = content.get_id()
            if isinstance(cid, str) and re.match('^\w*?\d+$', cid) is None:
                new_id = 'FILE%06d' % i
                i += 1
                self.id_maps[new_id] = cid
                content.change_id(new_id)
                new_fpath = content.copy_to(self.input_dir)
            else:
                # Put the content in the appropriate directory.
                new_fpath = content.copy_to(self.input_dir)
            self.num_input += 1
            logger.debug('%s saved for reading by reach.'
                         % new_fpath)
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
        for prefix in json_prefixes:
            base_prefix = path.basename(prefix)
            if base_prefix.isdecimal():
                base_prefix = int(base_prefix)
            elif base_prefix in self.id_maps.keys():
                base_prefix = self.id_maps[base_prefix]
            try:
                content = self._join_json_files(prefix, clear=True)
            except Exception as e:
                logger.exception(e)
                logger.error("Could not load result for prefix %s." % prefix)
                content = None
            self.add_result(base_prefix, content)
            logger.debug('Joined files for prefix %s.' % base_prefix)
        return self.results

    def clear_input(self):
        """Remove all the input files (at the end of a reading)."""
        for item in listdir(self.input_dir):
            item_path = path.join(self.input_dir, item)
            if path.isfile(item_path):
                remove(item_path)
                logger.debug('Removed input %s.' % item_path)
        return

    def _read(self, read_list, verbose=False, log=False):
        """Read the content, returning a list of ReadingData objects."""
        ret = []
        mem_tot = get_mem_total()
        if mem_tot is not None and mem_tot <= self.REACH_MEM + self.MEM_BUFFER:
            logger.error(
                "Too little memory to run reach. At least %s required." %
                (self.REACH_MEM + self.MEM_BUFFER)
            )
            logger.info("REACH not run.")
            return ret

        # Prep the content
        self.prep_input(read_list)

        # Make sure there is something to read before we start up Reach.
        if not self.num_input:
            return ret

        # Run REACH!
        logger.info("Beginning reach.")
        args = [
            'java',
            '-Dconfig.file=%s' % self.conf_file_path,
            '-jar', self.exec_path
        ]
        p = subprocess.Popen(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        # Monitor the logs and wait for Reach to finish.
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

        # Get the output
        ret = self.get_output()
        self.clear_input()

        return ret

    @staticmethod
    def get_processor(content):
        json_str = json.dumps(content)
        return reach.process_json_str(json_str)
