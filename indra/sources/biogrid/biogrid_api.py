import rdflib
import logging
from indra.sources.biogrid.processor import BiogridProcessor

logger = logging.getLogger('cwms_rdf')


def process_file(filename):
    """Processes a biogrid tab-separated filer.

    Parameters:
    -----------
    text: str
        The filename of the biogrid file

    Returns
    -------
    bp: indra.sources.cwms.BiogridProcessor
        A BiogridProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    bp = BiogridProcessor(filename)
    return bp
    
if __name__ == '__main__':
    fname = '/Users/daniel/Downloads/BIOGRID/BIOGRID-ALL-3.4.158.tab2.txt';bp = process_file(fname)
