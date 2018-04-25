from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import tempfile
import shutil
from indra.sources.isi.processor import IsiProcessor
import subprocess
import logging

logger = logging.getLogger('isi')

def process_directory(dir_name):
    """Process a directory filled with text files (it is sufficient for them
    to have extension .txt and be in the top level.
    
    Parameters
    ----------
    dir_name: str
        The directory containing input files from which to extract statements

    Returns
    -------
    ip: indra.sources.isi.processor.IsiProcessor
        A processor containing extracted statements
    """

    # Create a temporary directory to store the output
    output_dir = tempfile.mkdtemp('indra_isi_processor_output')
    tmp_dir = tempfile.mkdtemp('indra_isi_processor_tmp')

    # Form the command to invoke the ISI reader via Docker
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
    ip = IsiProcessor(output_dir)

    # Remove the temporary output directory
    shutil.rmtree(output_dir)
    shutil.rmtree(tmp_dir)

    return ip
