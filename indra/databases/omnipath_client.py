from __future__ import unicode_literals
from builtins import dict, str
import logging
import requests
from indra.databases import hgnc_client, uniprot_client
from indra.statements import *

logger = logging.getLogger("omnipath")

op_url = 'http://omnipathdb.org'

def _agent_from_up_id(up_id):
    """Build an Agent object from a Uniprot ID. Adds db_refs for both Uniprot
    and HGNC where available."""
    db_refs = {'UP': up_id}
    gene_name = uniprot_client.get_gene_name(up_id)
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
    ag = Agent(gene_name, db_refs=db_refs)
    return ag


def _stmts_from_op_mods(mod_list):
    """Build INDRA Statements from a list of Omnipath PTM list entries."""
    stmt_list = []
    for mod_entry in mod_list:
        enz = _agent_from_up_id(mod_entry['enzyme'])
        sub = _agent_from_up_id(mod_entry['substrate'])
        res = mod_entry['residue_type']
        pos = mod_entry['residue_offset']
        ref_list = mod_entry['references']
        evidence = [Evidence('omnipath', None, pmid)
                    for pmid in mod_entry['references']]
        if mod_entry['modification'] == 'phosphorylation':
            stmt = Phosphorylation(enz, sub, res, pos, evidence)
        stmt_list.append(stmt)
    return stmt_list


def get_all_modifications():
    """Get all PTMs from Omnipath as INDRA Statements.

    Returns
    -------
    list of Statements
    """
    params = {'format': 'json', 'fields':['sources', 'references']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    return _stmts_from_op_mods(res.json())


def get_modifications(up_list):
    """Get all PTMs from Omnipath for a list of proteins

    Parameters
    ----------
    up_list : list
        A list of Uniprot IDs.

    Returns
    -------
    list of Statements
    """
    params = {'format': 'json', 'fields':['sources', 'references']}
    gene_str = ','.join(up_list)
    ptm_url = '%s/ptms/%s' % (op_url, gene_str)
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    return _stmts_from_op_mods(res.json())

