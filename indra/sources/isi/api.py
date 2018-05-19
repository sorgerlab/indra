from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import subprocess
import logging
import os
import tempfile
import shutil
from indra.sources.isi.processor import IsiProcessor
from indra.sources.isi.preprocessor import IsiPreprocessor

logger = logging.getLogger('isi')


def process_text(text, pmid=None):
    """Process a string using the ISI reader and extract INDRA statements.

    Parameters
    ----------
    text : str
        A string of biomedical information to process
    pmid : Optional[str]
        The pmid associated with this text (or None if not specified)

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
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


def process_nxml(nxml_filename, pmid=None, extra_annotations=None):
    """Process an NXML file using the ISI reader

    First converts NXML to plain text and preprocesses it, then runs the ISI
    reader, and processes the output to extract INDRA Statements.

    Parameters
    ----------
    nxml_filename : str
        nxml file to process
    pmid : Optional[str]
        pmid of this nxml file, to be added to the Evidence object of the
        extracted INDRA statements
    extra_annotations : Optional[dict]
        Additional annotations to add to the Evidence object of all extracted
        INDRA statements. Extra annotations called 'interaction' are ignored
        since this is used by the processor to store the corresponding
        raw ISI output.

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
        A processor containing extracted Statements
    """
    if extra_annotations is None:
        extra_annotations = {}

    # Create a temporary directory to store the proprocessed input
    pp_dir = tempfile.mkdtemp('indra_isi_pp_output')

    pp = IsiPreprocessor(pp_dir)
    extra_annotations = {}
    pp.preprocess_nxml_file(nxml_filename, pmid, extra_annotations)

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
    isi_preprocessor : indra.sources.isi.preprocessor.IsiPreprocessor
        Preprocessor object that has already preprocessed the documents we
        want to read and process with the ISI reader
    num_processes : Optional[int]
        Number of processes to parallelize over
    specified_output_dir : Optional[str]
        The directory into which to put reader output; if omitted or None,
        uses a temporary directory.

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
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
    # We call realpath on all these paths so that any symbolic links
    # are generated out - this is needed on Mac
    input_binding = os.path.realpath(dir_name) + ':/input:ro'
    output_binding = os.path.realpath(output_dir) + ':/output:rw'
    tmp_binding = os.path.realpath(tmp_dir) + ':/temp:rw'
    command = ['docker', 'run', '-it', '--rm',
               '-v', input_binding, '-v', output_binding, '-v', tmp_binding,
               'sahilgar/bigmechisi', './myprocesspapers.sh',
               '-c', str(num_processes)]

    # Invoke the ISI reader
    logger.info('Running command:')
    logger.info(' '.join(command))
    ret = subprocess.call(command)
    if ret != 0:
        logger.error('Docker returned non-zero status code')

    # Process ISI output
    ip = IsiProcessor(output_dir, isi_preprocessor.pmids,
                      isi_preprocessor.extra_annotations)

    # Remove the temporary output directory
    if specified_output_dir is None:
        shutil.rmtree(output_dir)
    shutil.rmtree(tmp_dir)

    return ip
