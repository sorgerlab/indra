from __future__ import unicode_literals
import logging
import requests
from collections import Counter
import pypath.intera as pp_intera
from pypath import main as pypath_main, data_formats
from indra.databases import hgnc_client, uniprot_client
from indra.statements import modtype_to_modclass, Agent, Evidence, Complex

logger = logging.getLogger(__file__)

op_url = 'http://omnipathdb.org'
pa = pypath_main.PyPath()
urls = {'interactions': op_url + '/interactions',
        'ptms': op_url + '/ptms'}


def _get_text_refs(article_id_list):
    text_refs = {}
    for ref in article_id_list:
        name = ref['idtype'].upper()
        try:
            id = int(ref['value'])
        except ValueError:
            id = ref['value']
        text_refs[name] = id
    return text_refs


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


def _complex_agents_from_op_complex(up_id_string):
    """up_ids is a string of the format COMPLEX:UP1_UP2_..."""
    # Return list of contained agents
    if 'COMPLEX' in up_id_string:
        if ' ' in up_id_string:
            up_id_string = up_id_string.split()[-1]
        return [_agent_from_up_id(up_id) for up_id in
                up_id_string.split('COMPLEX:')[1].split('_')]
    else:
        return [_agent_from_up_id(up_id_string)]


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


def _stmts_from_op_pypath_graph(pp):
    """Build Complex statements from an igraph of ligand-receptor interactions

    Parameters
    ----------
    pp : pypath.main.PyPath
        An instance of a PyPath object containing the network
        representing ligand-receptor interactions
    """
    stmt_list = []
    for s, t in pp.graph.get_edgelist():
        edge_obj = pp.get_edge(s, t)

        # Get participating agents
        if isinstance(pp.vs[s]['name'], pp_intera.Complex):
            # Catch the odd pypath.intera.Complex objects
            src_string = str(pp.vs[s]['name'])
        else:
            src_string = pp.vs[s]['name']
        source_agents = _complex_agents_from_op_complex(src_string)

        if isinstance(pp.vs[t]['name'], pp_intera.Complex):
            # Catch the odd pypath.intera.Complex objects
            trg_string = str(pp.vs[t]['name'])
        else:
            trg_string = pp.vs[t]['name']
        target_agents = _complex_agents_from_op_complex(trg_string)

        # Assemble agent list
        agent_list = []
        for agent in [*source_agents, *target_agents]:
            if agent not in agent_list:
                agent_list.append(agent)

        # Get article IDs by support
        for ref_name, ref_set in edge_obj['refs_by_source'].items():
            for ref_obj in ref_set:
                # Check for PMID
                if ref_obj.pmid:
                    pmid = ref_obj.pmid
                else:
                    pmid = None
                try:
                    ref_info = ref_obj.info()
                    if ref_info.get('uids'):
                        uid = ref_info['uids'][0]
                        # Text refs
                        text_refs = _get_text_refs(ref_info[uid]['articleids'])
                        text_refs['nlmuniqueid'] = ref_info[uid]['nlmuniqueid']
                        text_refs['ISSN'] = ref_info[uid]['issn']
                        text_refs['ESSN'] = ref_info[uid]['essn']
                    else:
                        text_refs = None
                except TypeError as e:
                    logger.warning('Failed to load info')
                    logger.exception(e)
                    text_refs = None

                # If both pmid and text_refs is None, skip this Complex
                if pmid is None and text_refs is None:
                    continue

                # Get annotations
                annotations = {'omnipath_source': ref_name}
                if ref_name == 'Ramilowski2015' and\
                        edge_obj['ramilowski_sources']:
                    annotations['ramilowski_sources'] =\
                        edge_obj['ramilowski_sources']
                if edge_obj['cellphonedb_type']:
                    annotations['cellphonedb_type'] =\
                        edge_obj['cellphonedb_type']
                evidence = Evidence('omnipath', None, pmid,
                                    annotations=annotations,
                                    text_refs=text_refs)
                stmt_list.append(Complex(agent_list, evidence))

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


def get_all_rlint_pypath(reload_resources=False, force=False):
    """Get all receptor ligand interactions from the omnipath pypath module

    Returns
    -------
    stmts : list[indra.statements.Statement]
        A list of indra statements"""
    if reload_resources:
        # Todo wipe the cache (stored in ~/.pypath/cache) clean and
        #  re-download the resources. Warn the user that it takes a lot of
        #  time to download it
        pass
    pa.init_network(data_formats.ligand_receptor)
    return _stmts_from_op_pypath_graph(pa)
