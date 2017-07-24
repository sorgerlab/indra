from __future__ import unicode_literals
from builtins import dict, str
import logging
import requests
from collections import Counter
from indra.databases import hgnc_client, uniprot_client
from indra.statements import modtype_to_modclass, Agent, Evidence

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
    """Build Modification Statements from a list of Omnipath PTM entries."""
    stmt_list = []
    unhandled_mod_types = []
    for mod_entry in mod_list:
        enz = _agent_from_up_id(mod_entry['enzyme'])
        sub = _agent_from_up_id(mod_entry['substrate'])
        res = mod_entry['residue_type']
        pos = mod_entry['residue_offset']
        ref_list = mod_entry['references']
        evidence = [Evidence('omnipath', None, pmid)
                    for pmid in mod_entry['references']]
        mod_type = mod_entry['modification']
        modclass = modtype_to_modclass.get(mod_type)
        if modclass is None:
            if mod_type == 'cleavage':
                print(mod_entry)
                print()
            unhandled_mod_types.append(mod_type)
            continue
        else:
            stmt = modclass(enz, sub, res, pos, evidence)
        stmt_list.append(stmt)
    print(Counter(unhandled_mod_types))
    return stmt_list

#'cleavage',
#'proteolytic cleavage',

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

