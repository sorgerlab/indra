import requests
import json


def get_mesh_id(nlm_mesh):
    
    url_nlm2mesh = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=mesh&term=%s&retmode=JSON' % nlm_mesh
    r = requests.get(url_nlm2mesh)
    res = r.json()
    mesh_id = res['esearchresult']['idlist'][0]
    return mesh_id


def get_pcid(mesh_id):

    url_mesh2pcid = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=mesh&id=%s&db=pccompound&retmode=JSON' % mesh_id
    r = requests.get(url_mesh2pcid)
    res = r.json()
    pcid = res['linksets'][0]['linksetdbs'][0]['links'][0]
    return pcid


def get_chembl_id(nlm_mesh):

    mesh_id = get_mesh_id(nlm_mesh)
    pcid = get_pcid(mesh_id)
    url_mesh2pcid = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/JSON' % pcid
    r = requests.get(url_mesh2pcid)
    res = r.json()
    synonyms = res['InformationList']['Information'][0]['Synonym']
    chembl_id = [syn for syn in synonyms if 'CHEMBL' in syn and 'SCHEMBL' not in syn][0]
    return chembl_id
