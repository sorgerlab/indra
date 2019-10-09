import os
import json
import pybel
import rdflib
import logging
from functools import lru_cache
from .processor import PybelProcessor


logger = logging.getLogger(__name__)


def process_pybel_neighborhood(gene_names, network_file=None,
                               network_type='belscript', **kwargs):
    """Return PybelProcessor around neighborhood of given genes in a network.

    This function processes the given network file and filters the returned
    Statements to ones that contain genes in the given list.

    Parameters
    ----------
    network_file : Optional[str]
        Path to the network file to process. If not given, by default, the
        BEL Large Corpus is used.
    network_type : Optional[str]
        This function allows processing both BEL Script files and JSON files.
        This argument controls which type is assumed to be processed, and the
        value can be either 'belscript' or 'json'. Default: bel_script

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    if network_file is None:
        # Use large corpus as base network
        network_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    os.path.pardir, os.path.pardir,
                                    os.path.pardir, 'data', 'large_corpus.bel')
    if network_type == 'belscript':
        bp = process_belscript(network_file, **kwargs)
    elif network_type == 'json':
        bp = process_json_file(network_file)

    filtered_stmts = []
    for stmt in bp.statements:
        found = False
        for agent in stmt.agent_list():
            if agent is not None:
                if agent.name in gene_names:
                    found = True
        if found:
            filtered_stmts.append(stmt)

    bp.statements = filtered_stmts

    return bp


def process_belrdf(rdf_str, print_output=True):
    """Deprecated: Return a BelRdfProcessor for a BEL/RDF string.

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
    from rdflib.plugins.parsers.ntriples import ParseError
    from .rdf_processor import BelRdfProcessor
    logger.warning('The BEL/RDF format is deprecated and the results of '
                   'this function are not guaranteed to be correct. '
                   'Running this function requires rdflib==4.2.1, which is '
                   'older than the rdflib dependency installed by default.')
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


@lru_cache(maxsize=100)
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
    if bp.annot_manager.failures:
        logger.warning('missing %d annotation pairs',
                       sum(len(v)
                           for v in bp.annot_manager.failures.values()))
    return bp


def process_belscript(file_name, **kwargs):
    """Return a PybelProcessor by processing a BEL script file.

    Key word arguments are passed directly to pybel.from_path,
    for further information, see
    pybel.readthedocs.io/en/latest/io.html#pybel.from_path
    Some keyword arguments we use here differ from the defaults
    of PyBEL, namely we set `citation_clearing` to False
    and `no_identifier_validation` to True.

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
    if 'citation_clearing' not in kwargs:
        kwargs['citation_clearing'] = False
    if 'no_identifier_validation' not in kwargs:
        kwargs['no_identifier_validation'] = True
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


def process_cbn_jgif_file(file_name):
    """Return a PybelProcessor by processing a CBN JGIF JSON file.

    Parameters
    ----------
    file_name : str
        The path to a CBN JGIF JSON file.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    with open(file_name, 'r') as jgf:
        return process_pybel_graph(pybel.from_cbn_jgif(json.load(jgf)))
