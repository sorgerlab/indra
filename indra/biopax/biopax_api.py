from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.java_vm import autoclass, JavaException
from indra.biopax import pathway_commons_client as pcc
from indra.biopax.processor import BiopaxProcessor

def process_pc_neighborhood(gene_names, neighbor_limit=1):
    """Returns a BiopaxProcessor for a PathwayCommons neighborhood query.

    The neighborhood query finds the neighborhood around a set of source genes.

    http://www.pathwaycommons.org/pc2/#graph

    http://www.pathwaycommons.org/pc2/#graph_kind

    Parameters
    ----------
    gene_names : list
        A list of HGNC gene symbols to search the neighborhood of.
        Examples: ['BRAF'], ['BRAF', 'MAP2K1']
    neighbor_limit : Optional[int]
        The number of steps to limit the size of the neighborhood around
        the gene names being queried. Default: 1

    Returns
    -------
    bp : BiopaxProcessor
        A BiopaxProcessor containing the obtained BioPAX model in bp.model.

    Notes
    -----
    As returned from this function, bp.statements is empty and one
    has to call methods of bp (e.g. bp.get_complexes()) to
    populate bp.statements.
    """
    model = pcc.graph_query('neighborhood', gene_names,
                            neighbor_limit=neighbor_limit)
    if model is not None:
        return process_model(model)


def process_pc_pathsbetween(gene_names, neighbor_limit=1):
    """Returns a BiopaxProcessor for a PathwayCommons paths-between query.

    The paths-between query finds the paths between a set of genes. Here
    source gene names are given in a single list and all directions of paths
    between these genes are considered.

    http://www.pathwaycommons.org/pc2/#graph

    http://www.pathwaycommons.org/pc2/#graph_kind

    Parameters
    ----------
    gene_names : list
        A list of HGNC gene symbols to search for paths between.
        Examples: ['BRAF', 'MAP2K1']
    neighbor_limit : Optional[int]
        The number of steps to limit the length of the paths between
        the gene names being queried. Default: 1

    Returns
    -------
    bp : BiopaxProcessor
        A BiopaxProcessor containing the obtained BioPAX model in bp.model.

    Notes
    -----
    As returned from this function, bp.statements is empty and one
    has to call methods of bp (e.g. bp.get_complexes()) to
    populate bp.statements.
    """
    model = pcc.graph_query('pathsbetween', gene_names,
                             neighbor_limit=neighbor_limit)
    if model is not None:
        return process_model(model)


def process_pc_pathsfromto(source_genes, target_genes, neighbor_limit=1):
    """Returns a BiopaxProcessor for a PathwayCommons paths-from-to query.

    The paths-from-to query finds the paths from a set of source genes to
    a set of target genes.

    http://www.pathwaycommons.org/pc2/#graph

    http://www.pathwaycommons.org/pc2/#graph_kind

    Parameters
    ----------
    source_genes : list
        A list of HGNC gene symbols that are the sources of paths being
        searched for.
        Examples: ['BRAF', 'RAF1', 'ARAF']
    target_genes : list
        A list of HGNC gene symbols that are the targets of paths being
        searched for.
        Examples: ['MAP2K1', 'MAP2K2']
    neighbor_limit : Optional[int]
        The number of steps to limit the length of the paths
        between the source genes and target genes being queried. Default: 1

    Returns
    -------
    bp : BiopaxProcessor
        A BiopaxProcessor containing the obtained BioPAX model in bp.model.

    Notes
    -----
    As returned from this function, bp.statements is empty and one
    has to call methods of bp (e.g. bp.get_complexes()) to
    populate bp.statements.
    """
    model = pcc.graph_query('pathsfromto', source_genes,
                             target_genes, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_owl(owl_filename):
    """Returns a BiopaxProcessor for a BioPAX OWL file.

    Parameters
    ----------
    owl_filename : string
        The name of the OWL file to process.

    Returns
    -------
    bp : BiopaxProcessor
        A BiopaxProcessor containing the obtained BioPAX model in bp.model.

    Notes
    -----
    As returned from this function, bp.statements is empty and one
    has to call methods of bp (e.g. bp.get_complexes()) to
    populate bp.statements.
    """
    model = pcc.owl_to_model(owl_filename)
    return process_model(model)


def process_model(model):
    """Returns a BiopaxProcessor for a BioPAX model object.

    Parameters
    ----------
    model : org.biopax.paxtools.model.Model
        A BioPAX model object.

    Returns
    -------
    bp : BiopaxProcessor
        A BiopaxProcessor containing the obtained BioPAX model in bp.model.

    Notes
    -----
    As returned from this function, bp.statements is empty and one
    has to call methods of bp (e.g. bp.get_complexes()) to
    populate bp.statements.
    """
    bproc = BiopaxProcessor(model)
    # bproc.get_complexes()
    # bproc.get_phosphorylation()
    # bproc.print_statements()
    return bproc
