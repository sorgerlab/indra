# -*- coding: utf-8 -*-
"""
This module provides a simplified API for invoking the Turku Event Extraction
System (TEES) on text and extracting INDRA statement from TEES output.

See publication:
Jari Björne, Sofie Van Landeghem, Sampo Pyysalo, Tomoko Ohta, Filip Ginter,
Yves Van de Peer, Sofia Ananiadou and Tapio Salakoski, PubMed-Scale Event
Extraction for Post-Translational Modifications, Epigenetics and Protein
Structural Relations. Proceedings of BioNLP 2012, pages 82-90, 2012.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.sources.tees.processor import TEESProcessor
from indra import get_config

import os.path
import logging
import codecs
import tempfile
import shutil
import subprocess
import glob
import gzip
import tarfile
import re

from indra.sources.tees.parse_tees import tees_parse_networkx_to_dot
import networkx.algorithms.dag as dag

__all__ = ['run_on_text', 'process_text', 'extract_output']

logger = logging.getLogger(__name__)

# If TEES isn't specified, we will check to see if any of these directories
# contain all of the files in tees_installation_files; if so, we'll assume
# that it is a TEES installation.
tees_candidate_paths = ['../TEES', '~/TEES', '~/Downloads/TEES']
tees_installation_files = ['batch.py', 'classify.py', 'train.py',
                           'visualize.py']
tees_installation_dirs = ['Classifiers', 'Detectors', 'Evaluators', 'Core']


def process_text(text, pmid=None, python2_path=None):
    """Processes the specified plain text with TEES and converts output to
    supported INDRA statements. Check for the TEES installation is the
    TEES_PATH environment variable, and configuration file; if not found,
    checks candidate paths in tees_candidate_paths. Raises an exception if
    TEES cannot be found in any of these places.

    Parameters
    ----------
    text : str
        Plain text to process with TEES
    pmid : str
        The PMID from which the paper comes from, to be stored in the Evidence
        object of statements. Set to None if this is unspecified.
    python2_path : str
        TEES is only compatible with python 2. This processor invokes this
        external python 2 interpreter so that the processor can be run in
        either python 2 or python 3. If None, searches for an executible named
        python2 in the PATH environment variable.

    Returns
    -------
    tp : TEESProcessor
        A TEESProcessor object which contains a list of INDRA statements
        extracted from TEES extractions
    """
    # Try to locate python2 in one of the directories of the PATH environment
    # variable if it is not provided
    if python2_path is None:
        for path in os.environ["PATH"].split(os.pathsep):
            proposed_python2_path = os.path.join(path, 'python2.7')
            if os.path.isfile(proposed_python2_path):
                python2_path = proposed_python2_path
                print('Found python 2 interpreter at', python2_path)
                break
    if python2_path is None:
        raise Exception('Could not find python2 in the directories ' +
                        'listed in the PATH environment variable. ' +
                        'Need python2 to run TEES.')

    # Run TEES
    a1_text, a2_text, sentence_segmentations = run_on_text(text,
                                                                python2_path)

    # Run the TEES processor
    tp = TEESProcessor(a1_text, a2_text, sentence_segmentations, pmid)
    return tp

def run_on_text(text, python2_path):
    """Runs TEES on the given text in a temporary directory and returns a
    temporary directory with TEES output.

    The caller should delete this directory when done with it. This function
    runs TEES and produces TEES output files but does not process TEES output
    into INDRA statements.

    Parameters
    ----------
    text : str
        Text from which to extract relationships
    python2_path : str
        The path to the python 2 interpreter

    Returns
    -------
    output_dir : str
        Temporary directory with TEES output. The caller should delete this
        directgory when done with it.
    """
    tees_path = get_config('TEES_PATH')

    if tees_path is None:
        # If TEES directory is not specifies, see if any of the candidate paths
        # exist and contain all of the files expected for a TEES installation.
        for cpath in tees_candidate_paths:
            cpath = os.path.expanduser(cpath)
            if os.path.isdir(cpath):
                # Check to see if it has all of the expected files and
                # directories
                has_expected_files = True
                for f in tees_installation_files:
                    fpath = os.path.join(cpath, f)
                    present = os.path.isfile(fpath)
                    has_expected_files = has_expected_files and present

                has_expected_dirs = True
                for d in tees_installation_dirs:
                    dpath = os.path.join(cpath, d)
                    present = os.path.isdir(dpath)
                    has_expected_dirs = has_expected_dirs and present

                if has_expected_files and has_expected_dirs:
                    # We found a directory with all of the files and
                    # directories  we expected in a TEES installation - let's
                    # assume it's a TEES installation
                    tees_path = cpath
                    print('Found TEES installation at ' + cpath)
                    break

    # Make sure the provided TEES directory exists
    if not os.path.isdir(tees_path):
        raise Exception('Provided TEES directory does not exist.')

    # Make sure the classify.py script exists within this directory
    classify_path = 'classify.py'
    # if not os.path.isfile(classify_path):
    #    raise Exception('classify.py does not exist in provided TEES path.')

    # Create a temporary directory to tag the shared-task files
    tmp_dir = tempfile.mkdtemp(suffix='indra_tees_processor')

    pwd = os.path.abspath(os.getcwd())

    try:
        # Write text to a file in the temporary directory
        text_path = os.path.join(tmp_dir, 'text.txt')
        # Had some trouble with non-ascii characters. A possible TODO item in
        # the future is to look into resolving this, for now just ignoring
        # non-latin-1 characters
        with codecs.open(text_path, 'w', encoding='latin-1', errors='ignore') \
                as f:
            f.write(text)

        # Run TEES
        output_path = os.path.join(tmp_dir, 'output')
        model_path = os.path.join(tees_path, 'tees_data/models/GE11-test/')
        command = [python2_path, classify_path, '-m', model_path,
                   '-i', text_path,
                   '-o', output_path]
        try:
            pwd = os.path.abspath(os.getcwd())
            os.chdir(tees_path)  # Change to TEES directory
            # print('cwd is:', os.getcwd())
            # out = subprocess.check_output(command, stderr=subprocess.STDOUT)
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, cwd=tees_path)
            p.wait()
            (so, se) = p.communicate()
            print(so)
            print(se)
            os.chdir(pwd)  # Change back to previous directory
            # print('cwd is:', os.getcwd())
            # print(out.decode('utf-8'))

        except BaseException as e:
            # If there's an error, print it out and then propagate the
            # exception
            os.chdir(pwd)  # Change back to previous directory
            # print (e.output.decode('utf-8'))
            raise e

    except BaseException as e:
        # If there was an exception, delete the temporary directory and
        # pass on the exception
        shutil.rmtree(tmp_dir)
        raise e
    # Return the temporary directory with the TEES output
    output_tuple = extract_output(tmp_dir)
    shutil.rmtree(tmp_dir)
    return output_tuple


def extract_output(output_dir):
    """Extract the text of the a1, a2, and sentence segmentation files from the
    TEES output directory. These files are located within a compressed archive.

    Parameters
    ----------
    output_dir : str
        Directory containing the output of the TEES system

    Returns
    -------
    a1_text : str
        The text of the TEES a1 file (specifying the entities)
    a2_text : str
        The text of the TEES a2 file (specifying the event graph)
    sentence_segmentations : str
        The text of the XML file specifying the sentence segmentation
    """

    # Locate the file of sentences segmented by the TEES system, described
    # in a compressed xml document
    sentences_glob = os.path.join(output_dir, '*-preprocessed*.xml.gz')
    sentences_filename_candidates = glob.glob(sentences_glob)

    # Make sure there is exactly one such file
    if len(sentences_filename_candidates) != 1:
        m = 'Looking for exactly one file matching %s but found %d matches'
        raise Exception(m % (
            sentences_glob, len(sentences_filename_candidates)))
        return None, None, None

    # Read in the sentence segmentation XML
    sentence_segmentation_filename = sentences_filename_candidates[0]
    with gzip.GzipFile(sentences_filename_candidates[0], 'r') as f:
        sentence_segmentations = f.read().decode('utf-8')

    # Create a temporary directory to which to extract the a1 and a2 files from
    # the tarball
    tmp_dir = tempfile.mkdtemp(suffix='indra_tees_processor')

    try:
        # Make sure the tarfile with the extracted events is in shared task
        # format is in the output directory
        tarfile_glob = os.path.join(output_dir, '*-events.tar.gz')
        candidate_tarfiles = glob.glob(tarfile_glob)
        if len(candidate_tarfiles) != 1:
            raise Exception('Expected exactly one match for glob %s' %
                            tarfile_glob)
            return None, None, None

        # Decide what tar files to extract
        # (We're not blindly extracting all files because of the security
        # warning in the documentation for TarFile.extractall
        # In particular, we want to make sure that the filename doesn't
        # try to specify a relative or absolute path other than the current
        # directory by making sure the filename starts with an alphanumeric
        # character.
        # We're also only interested in files with the .a1 or .a2 extension
        tar_file = tarfile.open(candidate_tarfiles[0])
        a1_file = None
        a2_file = None
        extract_these = []
        for m in tar_file.getmembers():
            if re.match('[a-zA-Z0-9].*.a[12]', m.name):
                extract_these.append(m)

                if m.name.endswith('.a1'):
                    a1_file = m.name
                elif m.name.endswith('.a2'):
                    a2_file = m.name
                else:
                    assert(False)

        # There should be exactly two files that match these criteria
        if len(extract_these) != 2 or a1_file is None or a2_file is None:
            raise Exception('We thought there would be one .a1 and one .a2' +
                            ' file in the tarball, but we got %d files total' %
                            len(extract_these))
            return None, None, None

        # Extract the files that we decided to extract
        tar_file.extractall(path=tmp_dir, members=extract_these)

        # Read the text of the a1 (entities) file
        with codecs.open(os.path.join(tmp_dir, a1_file), 'r',
                         encoding='utf-8') as f:
            a1_text = f.read()

        # Read the text of the a2 (events) file
        with codecs.open(os.path.join(tmp_dir, a2_file), 'r',
                         encoding='utf-8') as f:
            a2_text = f.read()

        # Now that we're done, remove the temporary directory
        shutil.rmtree(tmp_dir)

        # Return the extracted text
        return a1_text, a2_text, sentence_segmentations
    except BaseException as e:
        # If there was an exception, delete the temporary directory and
        # pass on the exception
        print('Not removing temporary directory: ' + tmp_dir)
        shutil.rmtree(tmp_dir)
        raise e
        return None, None, None
