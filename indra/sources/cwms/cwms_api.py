import rdflib
import logging
from indra.sources.cwms.processor import CWMSProcessor

logger = logging.getLogger('cwms')

def process_rdf_file(text, rdf_filename):
    """Process CWMS's RDF output for the given statement and returns a
    processor populated with INDRA statements.

    Parameters
    ----------
    text : str
        Sentence to process
    rdf_filename : str
        The RDF filename to process

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    cp = CWMSProcessor(text, rdf_filename)
    return cp

