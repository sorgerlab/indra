from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import os
import json
import shutil
import logging
import tempfile
import subprocess
from indra.sources.isi.processor import IsiProcessor
from indra.sources.isi.preprocessor import IsiPreprocessor

logger = logging.getLogger('isi')


def process_text(text, pmid=None, cleanup=True):
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

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def process_nxml(nxml_filename, pmid=None, extra_annotations=None,
                 cleanup=True):
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

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def process_preprocessed(isi_preprocessor, num_processes=1,
                         output_dir=None, cleanup=True):
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
            extra_annotations=isi_preprocessor.extra_annotations.get(fname, {}))
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


def process_output_folder(folder_path, pmids=None, extra_annotations=None):
    """Recursively extracts statements from all ISI output files in the
    given directory and subdirectories.

    Parameters
    ----------
    folder_path : str
        The directory to traverse
    """
    ips = []
    for entry in os.listdir(folder_path):
        full_entry_path = os.path.join(folder_path, entry)

        if os.path.isdir(full_entry_path):
            logger.warning('ISI processor: did not expect any ' +
                           'subdirectories in the output directory.')
            ip = process_output_folder(full_entry_path)
            ips.append(ip)
        elif entry.endswith('.json'):
            # Extract the corresponding file id
            m = re.match('([0-9]+)\.json', entry)
            if m is None:
                logger.warning('ISI processor:', entry, ' does not ' +
                               ' match expected format for output files.')
                pmid = None
                extra_annotations = {}
            else:
                doc_id = int(m.group(1))
                pmid = pmids.get(doc_id)
                extra_annotations = extra_annotations.get(doc_id)
            ip = process_json_file(full_entry_path, pmid, extra_annotations)
            ips.append(ip)
        else:
            logger.warning('ISI processor: did not expect any non-json ' +
                           'files in the output directory')
    if len(ips) > 1:
        for ip in ips[1:]:
            ips[0].statements += ip.statements

    if ips:
        return ips[0]
    else:
        return None


def process_json_file(file_path, pmid=None, extra_annotations=None):
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
    """
    logger.info('Extracting from %s' % file_path)
    with open(file_path, 'rb') as fh:
        jd = json.load(fh)
        ip = IsiProcessor(jd, pmid, extra_annotations)
        ip.get_statements()
        return ip
