from collections import defaultdict

from .processor import *


logger = logging.getLogger(__name__)


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
    mp = process_directory(directory_name, lazy=True)

    for statement in mp.iter_statements():
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict


def process_directory(directory_name, lazy=False):
    """Processes a directory filled with CSXML files, first normalizing the
    character encodings to utf-8, and then processing into a list of INDRA
    statements.

    Parameters
    ----------
    directory_name : str
        The name of a directory filled with csxml files to process
    lazy : bool
        If True, the statements will not be generated immediately, but rather
        a generator will be formulated, and statements can be retrieved by
        using `iter_statements`. If False, the `statements` attribute will be
        populated immediately. Default is False.

    Returns
    -------
    mp : indra.sources.medscan.processor.MedscanProcessor
        A MedscanProcessor populated with INDRA statements extracted from the
        csxml files
    """

    # Parent Medscan processor containing extractions from all files
    mp = MedscanProcessor()
    mp.process_directory(directory_name, lazy)
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
    mp = process_file(file_name, lazy=True)

    for statement in mp.iter_statements():
        s_dict[statement.evidence[0].pmid].append(statement)
    return s_dict


def process_file(filename, interval=None, lazy=False):
    """Process a CSXML file for its relevant information.

    Consider running the fix_csxml_character_encoding.py script in
    indra/sources/medscan to fix any encoding issues in the input file before
    processing.

    Attributes
    ----------
    filename : str
        The csxml file, containing Medscan XML, to process
    interval : (start, end) or None
        Select the interval of documents to read, starting with the
        `start`th document and ending before the `end`th document. If
        either is None, the value is considered undefined. If the value
        exceeds the bounds of available documents, it will simply be
        ignored.
    lazy : bool
        If True, the statements will not be generated immediately, but rather
        a generator will be formulated, and statements can be retrieved by
        using `iter_statements`. If False, the `statements` attribute will be
        populated immediately. Default is False.

    Returns
    -------
    mp : MedscanProcessor
        A MedscanProcessor object containing extracted statements
    """
    mp = MedscanProcessor()
    mp.process_csxml_file(filename, interval, lazy)
    return mp
