import codecs
import logging
import lxml.etree
from collections import namedtuple, defaultdict
from indra.sources.medscan.processor import *
import os
import glob
import tempfile
import shutil
from math import floor
import pickle
from indra.sources.medscan.fix_csxml_character_encoding import \
        fix_character_encoding
import time

logger = logging.getLogger('medscan')

def process_directory_statements_sorted_by_pmid(directory_name,
                                                medscan_resource_dir):
    """Processes a directory filled with CSXML files. For each file, first
    normalizes the character encoding to utf-8, and then processes into
    INDRA statements. Returns a dictionary mapping pmids to a list of
    statements corresponding to that pmid."""
    s_dict = defaultdict(list)
    mp = process_directory(directory_name, medscan_resource_dir)

    for statement in mp.statements:
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict

def process_directory(directory_name, medscan_resource_dir):
    """Processes a directory filled with CSXML files. For each file, first
    normalizes the character encoding to utf-8, and then processes into
    INDRA statements."""

    # Parent Medscan processor containing extractions from all files
    mp = MedscanProcessor(medscan_resource_dir)
    mp.log_entities = defaultdict(int)

    # Create temporary directory into which to put the csxml files with
    # normalized character encodings
    tmp_dir = tempfile.mkdtemp('indra_medscan_processor')
    tmp_file = os.path.join(tmp_dir, 'fixed_char_encoding')

    # Process each file
    glob_pattern = os.path.join(directory_name, '*.csxml')
    files = glob.glob(glob_pattern)
    num_files = float(len(files))
    print(int(num_files), 'files to read')
    percent_done = 0
    files_processed = 0
    start_time_s = time.time()

    for filename in files:
        fix_character_encoding(filename, tmp_file)
        mp_file = process_file(tmp_file, medscan_resource_dir, None, False)

        mp.statements.extend(mp_file.statements)
        mp.num_entities += mp_file.num_entities
        mp.num_entities_not_found += mp_file.num_entities_not_found
        if mp.unmapped_urns is None:
            mp.unmapped_urns = mp_file.unmapped_urns
        else:
            mp.unmapped_urns = mp.unmapped_urns.update(mp_file.unmapped_urns)

        for k in mp_file.log_entities:
            mp.log_entities[k] = mp.log_entities[k] + mp_file.log_entities[k]

        percent_done_now = floor(100.0 * files_processed / num_files)
        if percent_done_now > percent_done:
            percent_done = percent_done_now
            ellapsed_s = time.time() - start_time_s
            ellapsed_min = ellapsed_s / 60.0

            msg = 'Processed %d of %d files (%f%% complete, %f minutes)' % \
                    (files_processed, num_files, percent_done, ellapsed_min)
            print(msg)
        files_processed += 1

    # Delete the temporary directory
    shutil.rmtree(tmp_dir)


    return mp

def process_file_sorted_by_pmid(file_name, medscan_resource_dir):
    """Processed a file and returns a dictionary mapping pmids to a list of
    statements corresponding to that pmid."""
    s_dict = defaultdict(list)
    mp = process_file(file_name, medscan_resource_dir)

    for statement in mp.statements:
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict


def process_file(filename, medscan_resource_dir, num_documents=None,
                 progress_updates=True):
    """Process a CSXML file for its relevant information.

    Consider running the fix_csxml_character_encoding.py script in
    indra/sources/medscan to fix any encoding issues in the input file before
    processing.

    The CSXML format consists of a top-level `<batch>` root element containing
    a series of `<doc>` (document) elements, in turn containing `<sec>`
    (section) elements, and in turn containing `<sent>` (sentence) elements.

    Within the `<sent>` element, a series of additional elements appear
    in the following order:

    * `<toks>`, which contains a tokenized form of the sentence in its
      text attribute
    * `<textmods>`, which describes any preprocessing/normalization done to
      the underlying text
    * `<match>` elements, each of which contains one of more `<entity>`
      elements, describing entities in the text with their identifiers.
      The local IDs of each entities are given in the `msid` attribute of
      this element; these IDs are then referenced in any subsequent SVO
      elements.
    * `<svo>` elements, representing subject-verb-object triples. SVO elements
      with a `type` attribute of `CONTROL` represent normalized regulation
      relationships; they often represent the normalized extraction of the
      immediately preceding (but unnormalized SVO element). However, in some
      cases there can be a "CONTROL" SVO element without its parent immediately
      preceding it.

    Attributes
    ----------
    filename : str
        The csxml file, containing Medscan XML, to process
    medscan_resource_dir : str
        A directory containing Unmapped Complexes.rnef and
        Unmapped Functional classes.rnef which describe unmapped URNs, or None
        if not available. These files are currently parsed but not used.
    num_documents : int
        The number of documents to process, or None to process all of the
        documents within the csxml file.
    progress_updates: bool
        Whether to print progress updates to the console

    Returns
    -------
    mp : MedscanProcessor
        A MedscanProcessor object containing extracted statements
    """
    mp = MedscanProcessor(medscan_resource_dir)

    if progress_updates:
        logger.info("Parsing %s to XML" % filename)
    with open(filename, 'rb') as f:
        mp.process_csxml_from_file_handle(f, num_documents)
    return mp
