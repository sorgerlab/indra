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


def process_directory(isi_preprocessor):
    """Process a directory filled with text files (it is sufficient for them
    to have extension .txt and be in the top level.

    Parameters
    ----------
    isi_preprocessor: indra.sources.isi.preprocessor.IsiPreprocessor
        Preprocessor object that has already preprocessed the documents we
        want to read and process with the ISI reader

    Returns
    -------
    ip: indra.sources.isi.processor.IsiProcessor
        A processor containing extracted statements
    """

    # Create a temporary directory to store the output
    output_dir = tempfile.mkdtemp('indra_isi_processor_output')
    tmp_dir = tempfile.mkdtemp('indra_isi_processor_tmp')

    # Form the command to invoke the ISI reader via Docker
    dir_name = isi_preprocessor.preprocessed_dir
    input_binding = dir_name + ':/input:ro'
    output_binding = output_dir + ':/output:rw'
    tmp_binding = tmp_dir + ':/temp:rw'
    command = ['docker', 'run', '-it', '--rm',
               '-v', input_binding, '-v', output_binding, '-v', tmp_binding,
               'sahilgar/bigmechisi', './myprocesspapers.sh']

    # Invoke the ISI reader
    print('Running command:')
    print(' '.join(command))
    ret = subprocess.call(command)
    if ret != 0:
        logger.error('Docker returned non-zero status code')

    # Process ISI output
    ip = IsiProcessor(output_dir, isi_preprocessor)

    # Remove the temporary output directory
    shutil.rmtree(output_dir)
    shutil.rmtree(tmp_dir)

    return ip

def test_api():
    preprocessed_dir = tempfile.mkdtemp('indra_isi_preprocessed')
    preprocessor = IsiPreprocessor(preprocessed_dir)

    input_dir = '/Users/daniel/workspace/isi/test_input'
    for filename in os.listdir(input_dir):
        path = os.path.join(input_dir, filename)
        preprocessor.preprocess_plain_text_file(path, 12, {'foo': 'bar'})

    ip = process_directory(preprocessor)

    #shutil.rmtree(preprocessed_dir)
    print(preprocessed_dir)
    return ip
