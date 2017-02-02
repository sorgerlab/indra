from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import requests
from indra.databases import chebi_client


def process_query(drug, target):
    chebi_id = drug.db_refs['CHEBI']
    mesh_id = drug.db_refs['MESH']
    if chebi_id:
        drug_chembl_id = chebi_client.get_chembl_id(chebi_id)
    elif mesh_id:
        drug_chembl_id = get_chembl_id(mesh_id)

    target_uid = target.db_refs['UP']
    target_chembl_id = get_target_chemblid(target_uid)

    json_dict = send_query(drug_chembl_id, target_chembl_id)
    return process_json(json_dict)


def process_json(json_dict):
    assays = []
    for assay in json_dict['activities']:
        p = Activity(assay)
        p.pmid(assay)
        assays.append(p)
    return assays


def get_target_chemblid(target_uid):
    url = 'https://www.ebi.ac.uk/chembl/api/data/target.json?' + \
          'target_components__accession=%s' % target_uid
    r = requests.get(url)
    r.raise_for_status()
    js = r.json()
    target_chemblid = js['targets'][1]['target_chembl_id']
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
    chembl_id = [syn for syn in synonyms if 'CHEMBL' in syn
                 and 'SCHEMBL' not in syn][0]
    return chembl_id


def send_query(drug_chemblid, target_chemblid):
    url = 'https://www.ebi.ac.uk/chembl/api/data/activity.json'
    params = {'molecule_chembl_id': drug_chemblid,
              'target_chembl_id': target_chemblid}
    r = requests.get(url, params=params)
    r.raise_for_status()
    js = r.json()
    return js


class Activity(object):
    def __init__(self, assay):
        self.description = assay['assay_description']
        self.metric = assay['standard_type'] + ': ' + \
                      str(assay['standard_value']) + \
                      str(assay['standard_units'])

    def pmid(self, assay):
        chembl_doc_id = str(assay['document_chembl_id'])
        url_pmid = 'https://www.ebi.ac.uk/chembl/api/data/document.json'
        params = {'document_chembl_id': chembl_doc_id}
        r = requests.get(url_pmid, params=params)
        r.raise_for_status()
        js = r.json()
        self.pmid = 'PMID' + str(js['documents'][0]['pubmed_id'])

    def __repr__(self):
        return '<%s reported in %s>' % (self.metric, self.pmid)