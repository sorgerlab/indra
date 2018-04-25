import logging
import nltk
import os

logger = logging.getLogger('isi')

class IsiPreprocessor(object):
    """Preprocesses a set of documents, one by one, and adds the preprocessed
    text to a temporary directory in a format suitable for the ISI reader.
    The ISI reader requires plain text with one sentence per line.

    Attributes
    ----------
    preprocessed_dir: str
        The directory holding the literature text preprocessed and sentence
        tokenized in a format suitable for the ISI reader
    next_file_id: int
        The next file with preprocessed text will be named next_file_id.txt
    pmids: dict
        A dictionary mapping file ids to the pmid of the text corresponding
        to that file, can be None if unknown
    extra_annotations: dict
        A dictionary mapping file ids to a (possibly empty) dictionary with
        additional annotations to include for statements extracted from this
        document
    """

    def __init__(self, preprocessed_dir):
        self.preprocessed_dir = preprocessed_dir
        self.next_file_id = 1
        self.pmids = {}
        self.extra_annotations = {}

        # This directory should be empty
        contents = os.listdir(preprocessed_dir)
        if len(contents) != 0:
            logger.warning('IsiPreprocessor should get an empty directory in' +
                           ' which to store preprocessed files.')

    def preprocess_plain_text_string(self, text, pmid, extra_annotations):
        """Preprocesses plain text by tokenizing into sentences and writing
        each sentence on its own line in a plain text file. All other
        preprocessing functions ultimately call this one."""
        output_file = str(self.next_file_id) + '.txt'
        output_file = os.path.join(self.preprocessed_dir, output_file)

        # Tokenize sentence
        sentences = nltk.sent_tokenize(text)

        # Write sentences to text file
        with open(output_file, 'w') as f:
            for sentence in sentences:
                f.write(sentence + '\n')

        # Store annotations
        self.pmids[self.next_file_id] = pmid
        self.extra_annotations[self.next_file_id] = extra_annotations
        
        # Increment file id
        self.next_file_id += 1

    def preprocess_plain_text_file(self, filename, pmid, extra_annotations):
        with open(filename, 'r') as f:
            content = f.read()
            self.preprocess_plain_text_string(content, pmid,
                                              extra_annotations)

    def preprocess_nxml_string(self, nxml, pmid, extra_annotations):
        raise NotImplementedError()

    def preprocess_nxml_file(self, filename, pmid, extra_annotations):
        raise NotImplementedError()
