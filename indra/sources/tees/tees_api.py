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

import os.path
import logging
import codecs

from indra.sources.tees.parse_tees import tees_parse_networkx_to_dot
import networkx.algorithms.dag as dag

logger = logging.getLogger('tees')

# If TEES isn't specified, we will check to see if any of these directories
# contain all of the files in tees_installation_files; if so, we'll assume
# that it is a TEES installation.
tees_candidate_paths = ['../TEES', '~/TEES', '~/Downloads/TEES']
tees_installation_files = ['batch.py', 'classify.py', 'train.py',
                           'visualize.py']
tees_installation_dirs = ['Classifiers', 'Detectors', 'Evaluators', 'Core']


def process_text(text, pmid=None, tees_path=None, python2_path=None):
    """Processes the specified plain text with TEES and converts output to
    supported INDRA statements.

    Parameters
    ----------
    text: str
        Plain text to process with TEES
    pmid: str
        The PMID from which the paper comes from, to be stored in the Evidence
        object of statements. Set to None if this is unspecified.
    tees_path: str
        The path of the TEES installation directory containing classify.py.
        If None, searches several common paths.
    python2_path: str
        TEES is only compatible with python 2. This processor invokes this
        external python 2 interpreter so that the processor can be run in
        either python 2 or python 3. If None, searches for an executible named
        python2 in the PATH environment variable.

    Returns
    -------
    tp: TEESProcessor
        A TEESProcessor object which contains a list of INDRA statements
        extracted from TEES extractions
    """

    # If TEES directory is not specifies, see if any of the candidate paths
    # exist and contain all of the files expected for a TEES installation.
    for cpath in tees_candidate_paths:
        cpath = os.path.expanduser(cpath)
        if os.path.isdir(cpath):
            # Check to see if it has all of the expected files and directories
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
                # We found a directory with all of the files and directories
                # we expected in a TEES installation - let's assume it's a
                # TEES installation
                tees_path = cpath
                print('Found TEES installation at ' + cpath)
                break

    # If tees_path is None then we didn't find any installations
    if tees_path is None:
        raise Exception('Could not find TEES installation')

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

    # Run the TEES processor
    tp = TEESProcessor(text, pmid, tees_path, python2_path)
    return tp
