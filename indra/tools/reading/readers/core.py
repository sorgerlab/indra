import json
import logging
import tempfile

from .util import get_dir, get_time_stamp, formats

logger = logging.getLogger(__name__)

# Set a character limit for reach reading
CONTENT_CHARACTER_LIMIT = 5e5
CONTENT_MAX_SPACE_RATIO = 0.5


class ReadingData(object):
    """Object to contain the data produced by a reading.

    Parameters
    ----------
    content_id : int or str
        A unique identifier of the text content that produced the reading,
        which can be mapped back to that content.
    reader_class : type
        The class of the reader, a child of
        `indra.tools.reading.readers.core.Reader`.
    reader_version : str
        A string identifying the version of the underlying nlp reader.
    reading_format : str
        The format of the reading result. Options are in indra.db.formats.
    reading : str or dict
        The content of the reading result. A string in the format given by
        `reading_format`.
    """

    def __init__(self, content_id, reader_class, reader_version,
                 reading_format, reading):
        self.content_id = content_id
        self.reader_class = reader_class
        self.reader_version = reader_version
        self.format = reading_format
        self.reading = reading
        self._statements = None
        return

    def __repr__(self):
        return self.__class__.__name__ + "(content_id=%s, reader_class=%s)" \
               % (self.content_id, self.reader_class.__name__)

    def get_statements(self, reprocess=False, add_metadata=False):
        """General method to create statements."""
        if self._statements is None or reprocess:

            # Handle the case that there is no content.
            if self.reading is None:
                self._statements = []
                return []

            # Map to the different processors.
            processor = self.reader_class.get_processor(self.reading)

            # Get the statements from the processor, if it was resolved.
            if processor is None:
                logger.error("Production of statements from %s failed for %s."
                             % (self.reader_class.name, self.content_id))
                stmts = []
            else:
                stmts = processor.statements

            # Add some metadata to the annotations
            if add_metadata:
                meta_info = {'READER': self.reader_class.name.upper(),
                             'CONTENT_ID': self.content_id}
                self._statements = []
                for stmt in stmts:
                    stmt.evidence.text_refs.update(meta_info)
                    self._statements.append(stmt)

        return self._statements[:]

    def to_json(self):
        return {'content_id': self.content_id,
                'reader_name': self.reader_class.name,
                'reader_version': self.reader_version,
                'reading_format': self.format,
                'reading': self.reading}

    @classmethod
    def from_json(cls, jd):
        jd['reader_class'] = get_reader_class(jd.pop('reader_name'))
        stored_version = jd['reader_version']
        current_version = jd['reader_class'].get_version()
        if stored_version != current_version:
            logger.debug("Current reader version does not match stored "
                         "version: %s (current) vs %s (stored)"
                         % (current_version, stored_version))
        return cls(**jd)


class Reader(object):
    """This abstract object defines and some general methods for readers."""
    name = NotImplemented

    def __init__(self, base_dir=None, n_proc=1, check_content=True,
                 input_character_limit=CONTENT_CHARACTER_LIMIT,
                 max_space_ratio=CONTENT_MAX_SPACE_RATIO,
                 ResultClass=ReadingData):
        if base_dir is None:
            base_dir = self.name.lower() + '_run'
        self.n_proc = n_proc
        self.base_dir = get_dir(base_dir)
        tmp_dir = tempfile.mkdtemp(
            prefix='%s_job_%s' % (self.name.lower(), get_time_stamp()),
            dir=self.base_dir
        )
        self.tmp_dir = tmp_dir
        self.input_dir = get_dir(tmp_dir, 'input')
        self.id_maps = {}
        self.do_content_check = check_content
        self.input_character_limit = input_character_limit
        self.max_space_ratio = max_space_ratio
        self.results = []
        self.ResultClass = ResultClass
        return

    def __repr__(self):
        return self.__class__.__name__ + '(\'%s\', n_proc=%d)' \
               % (self.name, self.n_proc)

    def reset(self):
        self.results = []
        self.id_maps = {}
        return

    def add_result(self, content_id, content, **kwargs):
        """"Add a result to the list of results."""
        result_object = self.ResultClass(content_id, self.__class__,
                                         self.get_version(), formats.JSON,
                                         content, **kwargs)
        self.results.append(result_object)
        return

    def _check_content(self, content_str):
        """Check if the content is likely to be successfully read."""
        if self.do_content_check:
            space_ratio = float(content_str.count(' '))/len(content_str)
            if space_ratio > self.max_space_ratio:
                return "space-ratio: %f > %f" % (space_ratio,
                                                 self.max_space_ratio)
            if len(content_str) > self.input_character_limit:
                return "too long: %d > %d" % (len(content_str),
                                              self.input_character_limit)
        return None

    @classmethod
    def get_version(cls):
        raise NotImplementedError()

    def read(self, read_list, verbose=False, log=False):
        "Read a list of items and return a dict of output files."
        raise NotImplementedError()

    @staticmethod
    def get_processor(content):
        raise NotImplementedError()


class EmptyReader(Reader):
    """A class name to use for Readers that are not implemented yet."""


class ReadingError(Exception):
    pass


def get_reader_classes(parent=Reader):
    """Get all childless the descendants of a parent class, recursively."""
    children = parent.__subclasses__()
    descendants = children[:]
    for child in children:
        grandchildren = get_reader_classes(child)
        if grandchildren:
            descendants.remove(child)
            descendants.extend(grandchildren)
    return descendants


def get_reader_class(reader_name):
    """Get a particular reader class by name."""
    for reader_class in get_reader_classes():
        if reader_class.name.lower() == reader_name.lower():
            return reader_class
    else:
        logger.error("No such reader: %s" % reader_name)
        return None


def get_reader(reader_name, *args, **kwargs):
    """Get an instantiated reader by name."""
    return get_reader_class(reader_name)(*args, **kwargs)


def dump_readings(readings, filepath):
    """Dump a list of ReadingData objects to a file as JSON."""
    json_list = []
    for rd in readings:
        json_list.append(rd.to_json())

    with open(filepath, 'w') as f:
        json.dump(json_list, f)

    return


def load_readings(filepath):
    """Load readings from the given filepath."""
    with open(filepath, 'r') as f:
        json_list = json.load(f)
    return [ReadingData.from_json(jd) for jd in json_list]
