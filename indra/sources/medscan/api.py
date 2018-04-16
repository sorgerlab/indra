import codecs
import logging
import lxml.etree
from collections import namedtuple
from indra.sources.medscan.processor import *

logger = logging.getLogger('medscan')


def process_file(filename, medscan_resource_dir, num_documents=None):
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

    Returns
    -------
    mp : MedscanProcessor
        A MedscanProcessor object containing extracted statements
    """
    mp = MedscanProcessor(medscan_resource_dir)
    logger.info("Parsing %s to XML" % filename)
    with open(filename, 'rb') as f:
        mp.process_csxml_from_file_handle(f, num_documents)
    return mp
