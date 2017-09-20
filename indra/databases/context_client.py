from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import copy
from indra.databases import cbio_client
# Python 2
try:
    basestring
# Python 3
except:
    basestring = str


def get_protein_expression(gene_names, cell_types):
    """Return the protein expression levels of genes in cell types.

    Parameters
    ----------
    gene_names : list
        HGNC gene symbols for which expression levels are queried.
    cell_types : list
        List of cell type names in which expression levels are queried.
        The cell type names follow the CCLE database conventions.

        Example: LOXIMVI_SKIN, BT20_BREAST

    Returns
    -------
    res : dict[dict[float]]
        A dictionary keyed by cell line, which contains another dictionary
        that is keyed by gene name, with estimated protein amounts as values.
    """
    A = 0.2438361
    B = 3.0957627
    mrna_amounts = cbio_client.get_ccle_mrna(gene_names, cell_types)
    protein_amounts = copy(mrna_amounts)
    for cell_type in cell_types:
        amounts = mrna_amounts.get(cell_type)
        if amounts is None:
            continue
        for gene_name, amount in amounts.items():
            if amount is not None:
                protein_amount = 10**(A * amount + B)
                protein_amounts[cell_type][gene_name] = protein_amount
    return protein_amounts

def get_mutations(gene_names, cell_types):
    """Return protein amino acid changes in given genes and cell types.

    Parameters
    ----------
    gene_names : list
        HGNC gene symbols for which mutations are queried.
    cell_types : list
        List of cell type names in which mutations are queried.
        The cell type names follow the CCLE database conventions.

        Example: LOXIMVI_SKIN, BT20_BREAST

    Returns
    -------
    res : dict[dict[list]]
        A dictionary keyed by cell line, which contains another dictionary
        that is keyed by gene name, with a list of amino acid substitutions
        as values.
    """
    mutations = cbio_client.get_ccle_mutations(gene_names, cell_types)
    return mutations
