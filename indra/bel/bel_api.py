import rdflib
from rdflib.plugins.parsers.ntriples import ParseError

import ndex_client
from processor import BelProcessor

def process_ndex_neighborhood(gene_names, rdf_out='bel_output.rdf'):
    """Return a BelProcessor for an NDEx network neighborhood.

    Parameters
    ----------
    gene_names : list
        A list of HGNC gene symbols to search the neighborhood of.
        Example: ['BRAF', 'MAP2K1']
    rdf_out : Optional[str]
        Name of the output file to save the RDF returned by the web service.
        This is useful for debugging purposes or to repeat the same query
        on an offline RDF file later. Default: bel_output.rdf

    Returns
    -------
    bp : BelProcessor
        A BelProcessor object which contains INDRA Statements in bp.statements.

    Notes
    -----
    This function calls process_belrdf to the returned RDF string from the
    webservice.
    """
    network_id = '9ea3c170-01ad-11e5-ac0f-000c29cb28fb'
    #url_suffix = '/bel2rdf/v1/network/%s/asBELRDF/query' % network_id
    url_suffix = '/network/%s/asBELRDF/query' % network_id
    params = {'searchString': ' '.join(gene_names)}
    rdf = ndex_client.send_request(url_suffix, params)
    if rdf is None:
        print 'No response for NDEx neighborhood query.'
        return None
    with open(rdf_out, 'wt') as fh:
        fh.write(rdf.encode('utf-8'))
    bp = process_belrdf(rdf)
    return bp


def process_belrdf(rdf_str):
    """Return a BelProcessor for a BEL/RDF string.

    Parameters
    ----------
    rdf_str : str
        A BEL/RDF string to be processed. This will usually come from reading
        a .rdf file.

    Returns
    -------
    bp : BelProcessor
        A BelProcessor object which contains INDRA Statements in bp.statements.

    Notes
    -----
    This function calls all the specific get_type_of_mechanism()
    functions of the newly constructed BelProcessor to extract
    INDRA Statements.
    """
    g = rdflib.Graph()
    try:
        g.parse(data=rdf_str, format='nt')
    except ParseError:
        print 'Could not parse rdf.'
        return None
    # Build INDRA statements from RDF
    bp = BelProcessor(g)
    bp.get_complexes()
    bp.get_activating_subs()
    bp.get_modifications()
    bp.get_dephosphorylations()
    bp.get_activating_mods()
    bp.get_composite_activating_mods()
    bp.get_activity_activity()

    # Print some output about the process
    bp.print_statement_coverage()
    print "\n--- Converted INDRA Statements -------------"
    bp.print_statements()
    return bp
