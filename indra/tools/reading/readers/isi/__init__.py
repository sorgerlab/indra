import json

from indra.tools.reading.readers.core import Reader
from indra.tools.reading.readers.util import get_dir

from indra.sources.isi.api import run_isi, get_isi_version
from indra.sources.isi.processor import IsiProcessor
from indra.sources.isi.preprocessor import IsiPreprocessor


class IsiReader(Reader):

    name = 'ISI'

    def __init__(self, *args, **kwargs):
        super(IsiReader, self).__init__(*args, **kwargs)

        # Define some extra directories
        self.nxml_dir = get_dir(self.tmp_dir, 'nxmls')
        self.isi_temp_dir = get_dir(self.tmp_dir, 'temp')
        self.output_dir = get_dir(self.tmp_dir, 'output')

        return

    def read(self, read_list, verbose=False, log=False, n_per_proc=None):
        # Create a preprocessor
        pp = IsiPreprocessor(self.input_dir)

        # Preprocess all the content.
        for content in read_list:
            if content.is_format('nxml'):
                content.copy_to(self.nxml_dir)
                pp.preprocess_nxml_file(content.get_filepath(),
                                        content.get_id(), {})
            elif content.is_format('txt', 'text'):
                pp.preprocess_plain_text_string(content.get_text(),
                                                content.get_id(), {})
            else:
                raise ValueError("Invalid/unrecognized format: %s"
                                 % content.get_format())

        # Run ISI
        run_isi(self.input_dir, self.output_dir, self.isi_temp_dir,
                self.n_proc)

        # Process the outputs
        for fname, cid, extra_annots in pp.iter_outputs(self.output_dir):
            with open(fname, 'rb') as f:
                content = json.load(f)
            self.add_result(cid, content)

        return self.results

    @classmethod
    def get_version(cls):
        return get_isi_version()

    @staticmethod
    def get_processor(content):
        processor = IsiProcessor(content)
        processor.get_statements()
        return processor
