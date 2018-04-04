from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import rdflib
import logging
from rdflib.plugins.parsers.ntriples import ParseError

from indra.databases import ndex_client
from .belrdf_processor import BelRdfProcessor
from .pybel_processor import PybelProcessor
import pybel

logger = logging.getLogger('bel')

ndex_bel2rdf = 'http://bel2rdf.bigmech.ndexbio.org'

def process_ndex_neighborhood(gene_names, network_id=None,
                              rdf_out='bel_output.rdf', print_output=True):
    """Return a BelRdfProcessor for an NDEx network neighborhood.

    Parameters
    ----------
    gene_names : list
        A list of HGNC gene symbols to search the neighborhood of.
        Example: ['BRAF', 'MAP2K1']
    network_id : Optional[str]
        The UUID of the network in NDEx. By default, the BEL Large Corpus
        network is used.
    rdf_out : Optional[str]
        Name of the output file to save the RDF returned by the web service.
        This is useful for debugging purposes or to repeat the same query
        on an offline RDF file later. Default: bel_output.rdf

    Returns
    -------
    bp : BelRdfProcessor
        A BelRdfProcessor object which contains INDRA Statements in bp.statements.

    Notes
    -----
    This function calls process_belrdf to the returned RDF string from the
    webservice.
    """
    if network_id is None:
        network_id = '9ea3c170-01ad-11e5-ac0f-000c29cb28fb'
    url = ndex_bel2rdf + '/network/%s/asBELRDF/query' % network_id
    params = {'searchString': ' '.join(gene_names)}
    # The ndex_client returns the rdf as the content of a json dict
    res_json = ndex_client.send_request(url, params, is_json=True)
    if not res_json:
        logger.error('No response for NDEx neighborhood query.')
        return None
    if res_json.get('error'):
        error_msg = res_json.get('message')
        logger.error('BEL/RDF response contains error: %s' % error_msg)
        return None
    rdf = res_json.get('content')
    if not rdf:
        logger.error('BEL/RDF response is empty.')
        return None

    with open(rdf_out, 'wb') as fh:
        fh.write(rdf.encode('utf-8'))
    bp = process_belrdf(rdf, print_output=print_output)
    return bp


def process_belrdf(rdf_str, print_output=True):
    """Return a BelRdfProcessor for a BEL/RDF string.

    Parameters
    ----------
    rdf_str : str
        A BEL/RDF string to be processed. This will usually come from reading
        a .rdf file.

    Returns
    -------
    bp : BelRdfProcessor
        A BelRdfProcessor object which contains INDRA Statements in
        bp.statements.

    Notes
    -----
    This function calls all the specific get_type_of_mechanism()
    functions of the newly constructed BelRdfProcessor to extract
    INDRA Statements.
    """
    g = rdflib.Graph()
    try:
        g.parse(data=rdf_str, format='nt')
    except ParseError as e:
        logger.error('Could not parse rdf: %s' % e)
        return None
    # Build INDRA statements from RDF
    bp = BelRdfProcessor(g)
    bp.get_complexes()
    bp.get_activating_subs()
    bp.get_modifications()
    bp.get_activating_mods()
    bp.get_transcription()
    bp.get_activation()
    bp.get_conversions()

    # Print some output about the process
    if print_output:
        bp.print_statement_coverage()
        bp.print_statements()
    return bp


def process_pybel_graph(graph):
    """Return a PybelProcessor by processing a PyBEL graph.

    Parameters
    ----------
    graph : pybel.struct.BELGraph
        A PyBEL graph to process

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    bp = PybelProcessor(graph)
    bp.get_statements()
    return bp


def process_belscript(file_name, **kwargs):
    """Return a PybelProcessor by processing a BEL script file.

    Key word arguments are passed directly to pybel.from_path,
    for further information, see
    pybel.readthedocs.io/en/latest/io.html#pybel.from_path

    Parameters
    ----------
    file_name : str
        The path to a BEL script file.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    pybel_graph = pybel.from_path(file_name, **kwargs)
    return process_pybel_graph(pybel_graph)


def process_json_file(file_name):
    """Return a PybelProcessor by processing a Node-Link JSON file.

    For more information on this format, see:
    http://pybel.readthedocs.io/en/latest/io.html#node-link-json

    Parameters
    ----------
    file_name : str
        The path to a Node-Link JSON file.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    with open(file_name, 'rt') as fh:
        pybel_graph = pybel.from_json_file(fh, False)
    return process_pybel_graph(pybel_graph)
