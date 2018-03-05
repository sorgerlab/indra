import rdflib
import logging
from indra.sources.bbn import processor


logger = logging.getLogger('bbn')


def process_json_file(fname):
    """Process a JSON-LD file to extract Statements and return a processor.

    Parameters
    ----------
    fname : str
        The path to the JSON-LD file to be processed.

    Returns
    -------
    bp : indra.sources.bbn.BBNProcessor
        A BBNProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    graph = _load_graph(fname)
    bp = processor.BBNProcessor(graph)
    bp.get_statements()
    return bp


def _load_graph(fname):
    g = rdflib.Graph()
    with open(fname, 'rb') as fh:
        logger.info('Started loading graph from %s' % fname)
        g.parse(fh, format='json-ld')
        logger.info('Finished loading graph')
    return g
