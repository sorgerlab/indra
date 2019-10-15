from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import os
import glob
import json
import shutil
import logging
import tempfile
import subprocess
from indra.sources.isi.processor import IsiProcessor
from indra.sources.isi.preprocessor import IsiPreprocessor

logger = logging.getLogger(__name__)


def process_text(text, pmid=None, cleanup=True, add_grounding=True):
    """Process a string using the ISI reader and extract INDRA statements.

    Parameters
    ----------
    text : str
        A text string to process
    pmid : Optional[str]
        The PMID associated with this text (or None if not specified)
    cleanup : Optional[bool]
        If True, the temporary folders created for preprocessed reading input
        and output are removed. Default: True
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped

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
    ip = process_preprocessed(pp, add_grounding=add_grounding)

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def process_nxml(nxml_filename, pmid=None, extra_annotations=None,
                 cleanup=True, add_grounding=True):
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
    cleanup : Optional[bool]
        If True, the temporary folders created for preprocessed reading input
        and output are removed. Default: True
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped

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
    ip = process_preprocessed(pp, add_grounding=add_grounding)

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def process_preprocessed(isi_preprocessor, num_processes=1,
                         output_dir=None, cleanup=True, add_grounding=True):
    """Process a directory of abstracts and/or papers preprocessed using the
    specified IsiPreprocessor, to produce a list of extracted INDRA statements.

    Parameters
    ----------
    isi_preprocessor : indra.sources.isi.preprocessor.IsiPreprocessor
        Preprocessor object that has already preprocessed the documents we
        want to read and process with the ISI reader
    num_processes : Optional[int]
        Number of processes to parallelize over
    output_dir : Optional[str]
        The directory into which to put reader output; if omitted or None,
        uses a temporary directory.
    cleanup : Optional[bool]
        If True, the temporary folders created for preprocessed reading input
        and output are removed. Default: True
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
        A processor containing extracted statements
    """

    # Create a temporary directory to store the output
    if output_dir is None:
        output_dir = tempfile.mkdtemp('indra_isi_processor_output')
    else:
        output_dir = os.path.abspath(output_dir)
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

    ips = []
    for basename, pmid in isi_preprocessor.pmids.items():
        fname = os.path.join(output_dir, '%s.json' % basename)
        ip = process_json_file(fname, pmid=pmid,
            extra_annotations=isi_preprocessor.extra_annotations.get(fname, {}),
            add_grounding=add_grounding)
        ips.append(ip)

    # Remove the temporary output directory
    if output_dir is None:
        if cleanup:
            shutil.rmtree(output_dir)
        else:
            logger.info('Not cleaning up %s' % output_dir)
    if cleanup:
        shutil.rmtree(tmp_dir)
    else:
        logger.info('Not cleaning up %s' % output_dir)

    if len(ips) > 1:
        for ip in ips[1:]:
            ips[0].statements += ip.statements

    if ips:
        return ips[0]
    else:
        return None


def process_output_folder(folder_path, pmids=None, extra_annotations=None,
                          add_grounding=True):
    """Recursively extracts statements from all ISI output files in the
    given directory and subdirectories.

    Parameters
    ----------
    folder_path : str
        The directory to traverse
    pmids : Optional[str]
        PMID mapping to be added to the Evidence of the extracted INDRA
        Statements
    extra_annotations : Optional[dict]
        Additional annotations to add to the Evidence object of all extracted
        INDRA statements. Extra annotations called 'interaction' are ignored
        since this is used by the processor to store the corresponding
        raw ISI output.
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped
    """
    pmids = pmids if pmids is not None else {}
    extra_annotations = extra_annotations if \
        extra_annotations is not None else {}
    ips = []
    for entry in glob.glob(os.path.join(folder_path, '*.json')):
        entry_key = os.path.splitext(os.path.basename(entry))[0]
        # Extract the corresponding file id
        pmid = pmids.get(entry_key)
        extra_annotation = extra_annotations.get(entry_key)
        ip = process_json_file(entry, pmid, extra_annotation,
                               add_grounding=add_grounding)
        ips.append(ip)

    if len(ips) > 1:
        for ip in ips[1:]:
            ips[0].statements += ip.statements

    if ips:
        return ips[0]
    else:
        return None


def process_json_file(file_path, pmid=None, extra_annotations=None,
                      add_grounding=True):
    """Extracts statements from the given ISI output file.

    Parameters
    ----------
    file_path : str
        The ISI output file from which to extract statements
    pmid : int
        The PMID of the document being preprocessed, or None if not
        specified
    extra_annotations : dict
        Extra annotations to be added to each statement from this document
        (can be the empty dictionary)
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped
    """
    logger.info('Extracting from %s' % file_path)
    with open(file_path, 'rb') as fh:
        jd = json.load(fh)
        ip = IsiProcessor(jd, pmid, extra_annotations,
                          add_grounding=add_grounding)
        ip.get_statements()
        return ip
