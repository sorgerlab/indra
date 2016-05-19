import json
import logging
import requests

logger = logging.getLogger('context')

ndexbio_context = 'http://general.bigmech.ndexbio.org:8081/context/'

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
    res : str
        A json string containing the predicted protein expression levels of
        the given proteins in the given cell types as returned by the
        NDEx web service.
    """
    req_type = 'expression/cell_line'
    if isinstance(gene_names, basestring):
        gene_names = [gene_names]
    if isinstance(cell_types, basestring):
        cell_types = [cell_types]
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def get_mutations(gene_names, cell_types):
    """Return the mutation status of genes in cell types.

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
    res : str
        A json string containing the mutation status of
        the given proteins in the given cell types as returned by the
        NDEx web service.
    """
    req_type = 'mutation/cell_line'
    if isinstance(gene_names, basestring):
        gene_names = [gene_names]
    if isinstance(cell_types, basestring):
        cell_types = [cell_types]
    params = {g: cell_types for g in gene_names}
    res = send_request(req_type, params)
    return res

def send_request(req_type, params=None):
    """Send a request to the NDEx web service.

    Parameters
    ----------
    req_type : str
        The API endpoint for the NDEx web service.

        Example: expression/cell_line
    params : dict
        Dictionary of parameters as required by the NDEx API.

    Returns
    -------
    res_json : str
        A json string returned by the NDEx web service.
    """
    if params is None:
        params = {}
    res = requests.post(ndexbio_context + req_type, json=params)
    if res.status_code != 200:
        logger.error('Request to NDEx service returned with status %d' %
                     res.status_code)
        return None
    res_json = res.json()
    return res_json
