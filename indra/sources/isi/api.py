from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import tempfile
import shutil
from indra.sources.isi.processor import IsiProcessor
from indra.sources.isi.preprocessor import IsiPreprocessor
import subprocess
import logging
import os

logger = logging.getLogger('isi')

def process_text(text, pmid=None):
    """Processes a string using the ISI reader, performing processing and
    extracting.

    Parameters
    ----------
    text: str
        A string of biomedical information to process
    pmid: str
        The pmid associated with this text (or None if not specified)
    
    Returns
    -------
    ip: indra.sources.isi.processor.IsiProcessor
        A processor containing statements
    """
    # Create a temporary directory to store the proprocessed input
    pp_dir = tempfile.mkdtemp('indra_isi_pp_output')

    pp = IsiPreprocessor(pp_dir)
    extra_annotations = {}
    pp.preprocess_plain_text_string(text, pmid, extra_annotations)

    # Run the ISI reader and extract statements
    ip = process_preprocessed(pp)

    # Remove temporary directory with processed input
    shutil.rmtree(pp_dir)

    return ip

def process_preprocessed(isi_preprocessor, num_processes=1,
                         specified_output_dir=None):
    """Process a directory of abstracts and/or papers preprocessed using the
    specified IsiPreprocessor, to produce a list of extracted INDRA statements.

    Parameters
    ----------
    isi_preprocessor: indra.sources.isi.preprocessor.IsiPreprocessor
        Preprocessor object that has already preprocessed the documents we
        want to read and process with the ISI reader
    num_processes: int
        Number of processes to parallelize over
    specified_output_dir: str
        The directory into which to put reader output; if omitted or None,
        uses a temporary directory.

    Returns
    -------
    ip: indra.sources.isi.processor.IsiProcessor
        A processor containing extracted statements
    """

    # Create a temporary directory to store the output
    if specified_output_dir is None:
        output_dir = tempfile.mkdtemp('indra_isi_processor_output')
    else:
        output_dir = os.path.abspath(specified_output_dir)
    tmp_dir = tempfile.mkdtemp('indra_isi_processor_tmp')

    # Form the command to invoke the ISI reader via Docker
    dir_name = isi_preprocessor.preprocessed_dir
    input_binding = dir_name + ':/input:ro'
    output_binding = output_dir + ':/output:rw'
    tmp_binding = tmp_dir + ':/temp:rw'
    command = ['docker', 'run', '-it', '--rm',
               '-v', input_binding, '-v', output_binding, '-v', tmp_binding,
               'sahilgar/bigmechisi', './myprocesspapers.sh',
               '-c', str(num_processes)]

    # Invoke the ISI reader
    print('Running command:')
    print(' '.join(command))
    ret = subprocess.call(command)
    if ret != 0:
        logger.error('Docker returned non-zero status code')

    # Process ISI output
    ip = IsiProcessor(output_dir, isi_preprocessor)

    # Remove the temporary output directory
    if specified_output_dir is None:
        shutil.rmtree(output_dir)
    shutil.rmtree(tmp_dir)

    return ip


def test_api():
    preprocessed_dir = tempfile.mkdtemp('indra_isi_preprocessed')
    preprocessor = IsiPreprocessor(preprocessed_dir)

    input_dir = '/Users/daniel/workspace/isi/test_input'
    for filename in os.listdir(input_dir):
        path = os.path.join(input_dir, filename)
        preprocessor.preprocess_plain_text_file(path, 12,
                                                {'foo': 'bar',
                                                 'original_file': filename})

    ip = process_directory(preprocessor)

    # shutil.rmtree(preprocessed_dir)
    print(preprocessed_dir)
    return ip
