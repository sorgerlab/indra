from datetime import datetime

__all__ = ['process_text', 'process_nxml', 'process_preprocessed',
           'process_json_file', 'process_output_folder']

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


DOCKER_IMAGE_NAME = 'sahilgar/bigmechisi'
IN_ISI_DOCKER = os.environ.get('IN_ISI_DOCKER', 'false').lower() == 'true'


class IsiRuntimeError(Exception):
    pass


def process_text(text, pmid=None, **kwargs):
    """Process a string using the ISI reader and extract INDRA statements.

    Parameters
    ----------
    text : str
        A text string to process
    pmid : Optional[str]
        The PMID associated with this text (or None if not specified)
    num_processes : Optional[int]
        Number of processes to parallelize over
    cleanup : Optional[bool]
        If True, the temporary folders created for preprocessed reading input
        and output are removed. Default: True
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped
    molecular_complexes_only : Optional[bool]
        If True, only Complex statements between molecular entities are retained
        after grounding.

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
        A processor containing statements
    """
    cleanup = kwargs.get('cleanup', True)

    # Create a temporary directory to store the proprocessed input
    pp_dir = tempfile.mkdtemp('indra_isi_pp_output')

    pp = IsiPreprocessor(pp_dir)
    extra_annotations = {}
    pp.preprocess_plain_text_string(text, pmid, extra_annotations)

    # Run the ISI reader and extract statements
    ip = process_preprocessed(pp, **kwargs)

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def process_nxml(nxml_filename, pmid=None, extra_annotations=None, **kwargs):
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
    num_processes : Optional[int]
        Number of processes to parallelize over
    cleanup : Optional[bool]
        If True, the temporary folders created for preprocessed reading input
        and output are removed. Default: True
    add_grounding : Optional[bool]
        If True the extracted Statements' grounding is mapped
    molecular_complexes_only : Optional[bool]
        If True, only Complex statements between molecular entities are retained
        after grounding.

    Returns
    -------
    ip : indra.sources.isi.processor.IsiProcessor
        A processor containing extracted Statements
    """
    if extra_annotations is None:
        extra_annotations = {}

    cleanup = kwargs.get('cleanup', True)

    # Create a temporary directory to store the proprocessed input
    pp_dir = tempfile.mkdtemp('indra_isi_pp_output')

    pp = IsiPreprocessor(pp_dir)
    pp.preprocess_nxml_file(nxml_filename, pmid, extra_annotations)

    # Run the ISI reader and extract statements
    ip = process_preprocessed(pp, **kwargs)

    if cleanup:
        # Remove temporary directory with processed input
        shutil.rmtree(pp_dir)
    else:
        logger.info('Not cleaning up %s' % pp_dir)

    return ip


def _make_links(dirname, link_dir):
    """Make links to files in a directory.

    This is used when running from within the modified ISI docker.
    """
    # Create a directory in the root dir with the appropriate name.
    if not os.path.exists(link_dir):
        os.mkdir(link_dir)

    # Link each of the files in the directory to this new link_dir.
    for fname in os.listdir(dirname):
        if fname.startswith('.'):
            continue
        link = os.path.join(link_dir, fname)
        true = os.path.join(dirname, fname)
        os.symlink(true, link)

    return


def run_isi(input_dir, output_dir, tmp_dir, num_processes=1,
            verbose=True, log=False):
    base_command = ['/root/myprocesspapers.sh', '-c', str(num_processes)]

    if IN_ISI_DOCKER:
        _make_links(input_dir, '/input')
        os.mkdir('/output')
        os.mkdir('/temp')
        command = base_command
    else:
        # We call realpath on all these paths so that any symbolic links
        # are generated out - this is needed on Mac
        input_binding = os.path.realpath(input_dir) + ':/input:ro'
        output_binding = os.path.realpath(output_dir) + ':/output:rw'
        tmp_binding = os.path.realpath(tmp_dir) + ':/temp:rw'
        command = ['docker', 'run', '-it', '--rm',
                   '-v', input_binding, '-v', output_binding, '-v', tmp_binding,
                   DOCKER_IMAGE_NAME] + base_command

    # Invoke the ISI reader
    logger.info('Running command from within the docker:' if IN_ISI_DOCKER
                else 'Running command using the docker:')
    logger.info(' '.join(command))

    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)

    # Monitor the logs and wait for reading to end.
    log_file_str = ''
    for line in iter(p.stdout.readline, b''):
        log_line = 'ISI: ' + line.strip().decode('utf8')
        if verbose:
            logger.info(log_line)
        if log:
            log_file_str += log_line + '\n'

    if log:
        with open('isi_run.log', 'ab') as f:
            f.write(log_file_str.encode('utf8'))

    p_out, p_err = p.communicate()
    if p.returncode:
        logger.error('Problem running ISI:')
        logger.error('Stdout & Stderr: %s' % p_out.decode('utf-8'))
        raise IsiRuntimeError("Problem encountered running ISI.")

    logger.info("ISI finished.")

    if IN_ISI_DOCKER:
        _make_links('/output', output_dir)
        _make_links('/temp', tmp_dir)

    return


def get_isi_image_data():
    """Get the json data for the ISI docker image."""
    if IN_ISI_DOCKER:
        logger.error("Cannot read docker info from within the docker.")
        return {}

    ret = subprocess.run(['docker', 'image', 'inspect', DOCKER_IMAGE_NAME],
                         stdout=subprocess.PIPE)
    image_data = json.loads(ret.stdout)[0]
    return image_data


def get_isi_version():
    if IN_ISI_DOCKER:
        timestamp = os.path.getmtime('/root/myprocesspapers.sh')
        dt = datetime.fromtimestamp(timestamp)
    else:
        data = get_isi_image_data()
        dt = datetime.strptime(data['Created'].split('.')[0],
                               '%Y-%m-%dT%H:%M:%S')
    return dt.strftime('%Y%m%d')


def process_preprocessed(isi_preprocessor, num_processes=1,
                         output_dir=None, cleanup=True, add_grounding=True,
                         molecular_complexes_only=False):
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
    molecular_complexes_only : Optional[bool]
        If True, only Complex statements between molecular entities are retained
        after grounding.

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

    # Run the ISI reader
    run_isi(dir_name, output_dir, tmp_dir, num_processes)

    ips = []
    for fname, pmid, extra_annots in isi_preprocessor.iter_outputs(output_dir):
        ip = process_json_file(fname, pmid=pmid,
                               extra_annotations=extra_annots,
                               add_grounding=add_grounding,
                            molecular_complexes_only=molecular_complexes_only)
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
                          add_grounding=True, molecular_complexes_only=False):
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
    molecular_complexes_only : Optional[bool]
        If True, only Complex statements between molecular entities are retained
        after grounding.
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
                               add_grounding=add_grounding,
                            molecular_complexes_only=molecular_complexes_only)
        ips.append(ip)

    if len(ips) > 1:
        for ip in ips[1:]:
            ips[0].statements += ip.statements

    if ips:
        return ips[0]
    else:
        return None


def process_json_file(file_path, pmid=None, extra_annotations=None,
                      add_grounding=True, molecular_complexes_only=False):
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
    molecular_complexes_only : Optional[bool]
        If True, only Complex statements between molecular entities are retained
        after grounding.
    """
    logger.info('Extracting from %s' % file_path)
    with open(file_path, 'rb') as fh:
        jd = json.load(fh)
        ip = IsiProcessor(jd, pmid, extra_annotations,
                          add_grounding=add_grounding)
        ip.get_statements()
        if molecular_complexes_only:
            ip.retain_molecular_complexes()
        return ip
