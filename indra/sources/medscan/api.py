import os
import glob
import time
import shutil
import codecs
import pickle
import tempfile
import logging
import lxml.etree
from math import floor
from collections import namedtuple, defaultdict
from .processor import *
from .fix_csxml_character_encoding import fix_character_encoding


logger = logging.getLogger('medscan')


def process_directory_statements_sorted_by_pmid(directory_name):
    """Processes a directory filled with CSXML files, first normalizing the
    character encoding to utf-8, and then processing into INDRA statements
    sorted by pmid.

    Parameters
    ----------
    directory_name : str
        The name of a directory filled with csxml files to process

    Returns
    -------
    pmid_dict : dict
        A dictionary mapping pmids to a list of statements corresponding to
        that pmid
    """
    s_dict = defaultdict(list)
    mp = process_directory(directory_name)

    for statement in mp.statements:
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict


def process_directory(directory_name):
    """Processes a directory filled with CSXML files, first normalizing the
    character encodings to utf-8, and then processing into a list of INDRA
    statements.

    Parameters
    ----------
    directory_name : str
        The name of a directory filled with csxml files to process

    Returns
    -------
    mp : indra.sources.medscan.processor.MedscanProcessor
        A MedscanProcessor populated with INDRA statements extracted from the
        csxml files
    """

    # Parent Medscan processor containing extractions from all files
    mp = MedscanProcessor()
    mp.log_entities = defaultdict(int)

    # Create temporary directory into which to put the csxml files with
    # normalized character encodings
    tmp_dir = tempfile.mkdtemp('indra_medscan_processor')
    tmp_file = os.path.join(tmp_dir, 'fixed_char_encoding')

    # Process each file
    glob_pattern = os.path.join(directory_name, '*.csxml')
    files = glob.glob(glob_pattern)
    num_files = float(len(files))
    logger.info("%d files to read" % int(num_files))
    percent_done = 0
    files_processed = 0
    start_time_s = time.time()

    for filename in files:
        logger.info('Processing', filename)
        fix_character_encoding(filename, tmp_file)
        mp_file = process_file(tmp_file, None)

        mp.statements.extend(mp_file.statements)
        mp.num_entities += mp_file.num_entities
        mp.num_entities_not_found += mp_file.num_entities_not_found

        for k in mp_file.log_entities:
            mp.log_entities[k] = mp.log_entities[k] + mp_file.log_entities[k]

        percent_done_now = floor(100.0 * files_processed / num_files)
        if percent_done_now > percent_done:
            percent_done = percent_done_now
            ellapsed_s = time.time() - start_time_s
            ellapsed_min = ellapsed_s / 60.0

            msg = 'Processed %d of %d files (%f%% complete, %f minutes)' % \
                    (files_processed, num_files, percent_done, ellapsed_min)
            logger.info(msg)
        files_processed += 1

    # Delete the temporary directory
    shutil.rmtree(tmp_dir)

    return mp


def process_file_sorted_by_pmid(file_name):
    """Processes a file and returns a dictionary mapping pmids to a list of
    statements corresponding to that pmid.

    Parameters
    ----------
    file_name : str
        A csxml file to process

    Returns
    -------
    s_dict : dict
        Dictionary mapping pmids to a list of statements corresponding to
        that pmid
    """
    s_dict = defaultdict(list)
    mp = process_file(file_name)

    for statement in mp.statements:
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict


def process_file(filename, num_documents=None):
    """Process a CSXML file for its relevant information.

    Consider running the fix_csxml_character_encoding.py script in
    indra/sources/medscan to fix any encoding issues in the input file before
    processing.

    Attributes
    ----------
    filename : str
        The csxml file, containing Medscan XML, to process
    num_documents : int
        The number of documents to process, or None to process all of the
        documents within the csxml file.

    Returns
    -------
    mp : MedscanProcessor
        A MedscanProcessor object containing extracted statements
    """
    mp = MedscanProcessor()

    with open(filename, 'rb') as f:
        mp.process_csxml_from_file_handle(f, num_documents)
    return mp
