import os
import json
import pybel
import pickle
import logging
import requests
from functools import lru_cache
from .processor import PybelProcessor


logger = logging.getLogger(__name__)

large_corpus_url = ('https://github.com/cthoyt/selventa-knowledge/raw/master/'
                    'selventa_knowledge/large_corpus-20170611.bel.pickle')
small_corpus_url = ('https://github.com/cthoyt/selventa-knowledge/raw/master/'
                    'selventa_knowledge/'
                    'selventa-small-corpus-20150611.bel.pickle')


def process_small_corpus():
    """Return PybelProcessor with statements from Selventa Small Corpus.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        its statements attribute.
    """
    return process_pybel_network(network_type='graph_pickle_url',
                                 network_file=small_corpus_url)


def process_large_corpus():
    """Return PybelProcessor with statements from Selventa Large Corpus.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        its statements attribute.
    """
    return process_pybel_network(network_type='graph_pickle_url',
                                 network_file=large_corpus_url)


def process_pybel_network(network_type, network_file, **kwargs):
    """Return PybelProcessor by processing a given network file.

    Parameters
    ----------
    network_type : str
        The type of network that network_file is. The options are:
        belscript, json, cbn_jgif, graph_pickle, and graph_pickle_url.
    network_file : str
        Path to the network file/URL to process.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    if network_type == 'belscript':
        return process_belscript(network_file, **kwargs)
    elif network_type == 'json':
        return process_json_file(network_file)
    elif network_type == 'cbn_jgif':
        return process_cbn_jgif_file(network_file)
    elif network_type == 'graph_pickle_url':
        if not network_file:
            network_file = large_corpus_url
        res = requests.get(network_file)
        res.raise_for_status()
        graph = pickle.loads(res.content)
        return process_pybel_graph(graph)
    elif network_type == 'graph_pickle':
        with open(network_file, 'rb') as fh:
            graph = pickle.load(fh)
            return process_pybel_graph(graph)
    else:
        raise ValueError('Unknown network type: %s' % network_type)


def process_pybel_neighborhood(entity_names, network_type='graph_pickle_url',
                               network_file=None, **kwargs):
    """Return PybelProcessor around neighborhood of given genes in a network.

    This function processes the given network file and filters the returned
    Statements to ones that contain genes in the given list.

    Parameters
    ----------
    entity_names : list[str]
        A list of entity names (e.g., gene names) which will be used as the
        basis of filtering the result. If any of the Agents of an extracted
        INDRA Statement has a name appearing in this list, the Statement is
        retained in the result.
    network_type : Optional[str]
        The type of network that network_file is. The options are:
        belscript, json, cbn_jgif, graph_pickle, and graph_pickle_url.
        Default: graph_pickle_url
    network_file : Optional[str]
        Path to the network file/URL to process. If not given, by default, the
        Selventa Large Corpus is used via a URL pointing to a PyBEL Graph
        pickle.

    Returns
    -------
    bp : PybelProcessor
        A PybelProcessor object which contains INDRA Statements in
        bp.statements.
    """
    bp = process_pybel_network(network_type, network_file, **kwargs)
    filtered_stmts = []
    filter_names = set(entity_names)
    for stmt in bp.statements:
        found = False
        for agent in stmt.agent_list():
            if agent is not None:
                if agent.name in filter_names:
                    found = True
        if found:
            filtered_stmts.append(stmt)

    bp.statements = filtered_stmts

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


def process_belrdf(rdf_str, print_output=True):
    """Deprecated: Return a BelRdfProcessor for a BEL/RDF string.

    Parameters
    ----------
    rdf_str : str
        A BEL/RDF string to be processed. This will usually come from reading
        a .rdf file.
    print_output : Optional[bool]
        If True, print statistics of what has been extracted from the given
        BEL/RDF network. Default: True

    Returns
    -------
    bp : BelRdfProcessor
        A BelRdfProcessor object which contains INDRA Statements in
        its statements attribute.

    Notes
    -----
    This function calls all the specific get_type_of_mechanism()
    functions of the newly constructed BelRdfProcessor to extract
    INDRA Statements.
    """
    import rdflib
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


