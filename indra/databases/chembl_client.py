from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import logging
import requests
from sympy.physics import units
from indra.databases import chebi_client
from indra.statements import Inhibition, Agent, Evidence

logger = logging.getLogger('chembl')


def get_inhibition(drug, target):
    chebi_id = drug.db_refs.get('CHEBI')
    mesh_id = drug.db_refs.get('MESH')
    if chebi_id:
        drug_chembl_id = chebi_client.get_chembl_id(chebi_id)
    elif mesh_id:
        drug_chembl_id = get_chembl_id(mesh_id)
    else:
        logger.error('Drug missing ChEBI or MESH grounding.')
        return None

    target_upid = target.db_refs.get('UP')
    if not target_upid:
        logger.error('Target missing UniProt grounding.')
        return None
    target_chembl_id = get_target_chemblid(target_upid)

    logger.info('Drug: %s, Target: %s' % (drug_chembl_id, target_chembl_id))
    query_dict = {'query': 'activity',
                  'params': {'molecule_chembl_id': drug_chembl_id,
                             'target_chembl_id': target_chembl_id,
                             'limit': 10000}
                  }
    res = send_query(query_dict)
    evidence = []
    for assay in res['activities']:
        ev = get_evidence(assay)
        if not ev:
            continue
        evidence.append(ev)
    st = Inhibition(drug, target, evidence=evidence)
    return st


def send_query(query_dict):
    """Query ChEMBL API

    Parameters
    ----------
    query_dict : dict
        'query' : string of the endpoint to query
        'params' : dict of params for the query
    Returns
    -------
    js : dict
        dict parsed from json that is unique to the submitted query

    Example
    -------
    >>> query_dict = {'query': 'target',
    >>>               'params': {'target_chembl_id': 'CHEMBL5145',
    >>>               'limit': 1}}
    >>> send_query(query_dict)
    """
    query = query_dict['query']
    params = query_dict['params']
    url = 'https://www.ebi.ac.uk/chembl/api/data/' + query + '.json'
    r = requests.get(url, params=params)
    r.raise_for_status()
    js = r.json()
    return js


def query_target(target_chembl_id):
    """Query ChEMBL API target by id

    Parameters
    ----------
    target_chembl_id : str
    Returns
    -------
    target : dict
        dict parsed from json that is unique for the target
    """
    query_dict = {'query': 'target',
                  'params': {'target_chembl_id': target_chembl_id,
                             'limit': 1}}
    res = send_query(query_dict)
    assert(res['page_meta']['total_count'] == 1)
    target = res['targets'][0]
    return target


def activities_by_target(activities):
    """Get back lists of activities in a dict keyed by ChEMBL target id
    Parameters
    ----------
    activities : dict
        response from a query returning activities for a drug
    Returns
    -------
    targ_act_dict : dict
        dictionary keyed to ChEMBL target ids with lists of activity ids
    """
    targ_act_dict = defaultdict(lambda: [])
    for activity in activities:
        target_chembl_id = activity['target_chembl_id']
        activity_id = activity['activity_id']
        targ_act_dict[target_chembl_id].append(activity_id)
    for target_chembl_id in targ_act_dict:
        targ_act_dict[target_chembl_id] = \
            list(set(targ_act_dict[target_chembl_id]))
    return targ_act_dict


def get_protein_targets_only(target_chembl_ids):
    """Given list of ChEMBL target ids, return dict of only SINGLE PROTEIN targ
    Parameters
    ----------
    target_chembl_ids : list
        list of chembl_ids as strings
    Returns
    -------
    protein_targets : dict
        dictionary keyed to ChEMBL target ids with lists of activity ids
    """
    protein_targets = {}
    for target_chembl_id in target_chembl_ids:
        target = query_target(target_chembl_id)
        if 'SINGLE PROTEIN' in target['target_type']:
            protein_targets[target_chembl_id] = target
    return protein_targets


def get_evidence(assay):
    kin = get_kinetics(assay)
    source_id = assay.get('assay_chembl_id')
    if not kin:
        return None
    annotations = {'kinetics': kin}
    chembl_doc_id = str(assay.get('document_chembl_id'))
    pmid = get_pmid(chembl_doc_id)
    ev = Evidence(source_api='chembl', pmid=pmid, source_id=source_id,
                  annotations=annotations)
    return ev


def get_kinetics(assay):
    try:
        val = float(assay.get('standard_value'))
    except TypeError:
        logger.warning('Invalid assay value: %s' % assay.get('standard_value'))
        return None
    unit = assay.get('standard_units')
    if unit == 'nM':
        unit_sym = 1e-9 * units.mol / units.liter
    elif unit == 'uM':
        unit_sym = 1e-6 * units.mol / units.liter
    else:
        logger.warning('Unhandled unit: %s' % unit)
        return None
    param_type = assay.get('standard_type')
    if param_type not in ['IC50', 'Kd']:
        logger.warning('Unhandled parameter type: %s' % param_type)
        return None
    kin = {param_type: val * unit_sym}
    return kin


def get_pmid(doc_id):
    url_pmid = 'https://www.ebi.ac.uk/chembl/api/data/document.json'
    params = {'document_chembl_id': doc_id}
    res = requests.get(url_pmid, params=params)
    js = res.json()
    pmid = str(js['documents'][0]['pubmed_id'])
    return pmid


def get_target_chemblid(target_upid):
    url = 'https://www.ebi.ac.uk/chembl/api/data/target.json'
    params = {'target_components__accession': target_upid}
    r = requests.get(url, params=params)
    r.raise_for_status()
    js = r.json()
    target_chemblid = js['targets'][0]['target_chembl_id']
    return target_chemblid


def get_mesh_id(nlm_mesh):
    url_nlm2mesh = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {'db': 'mesh', 'term': nlm_mesh, 'retmode': 'JSON'}
    r = requests.get(url_nlm2mesh, params=params)
    res = r.json()
    mesh_id = res['esearchresult']['idlist'][0]
    return mesh_id


def get_pcid(mesh_id):
    url_mesh2pcid = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    params = {'dbfrom': 'mesh', 'id': mesh_id,
              'db': 'pccompound', 'retmode': 'JSON'}
    r = requests.get(url_mesh2pcid, params=params)
    res = r.json()
    pcid = res['linksets'][0]['linksetdbs'][0]['links'][0]
    return pcid


def get_chembl_id(nlm_mesh):
    mesh_id = get_mesh_id(nlm_mesh)
    pcid = get_pcid(mesh_id)
    url_mesh2pcid = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/' + \
                    'cid/%s/synonyms/JSON' % pcid
    r = requests.get(url_mesh2pcid)
    res = r.json()
    synonyms = res['InformationList']['Information'][0]['Synonym']
    chembl_id = [syn for syn in synonyms
                 if 'CHEMBL' in syn and 'SCHEMBL' not in syn][0]
    return chembl_id
