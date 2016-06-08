import chembl_client
from chembl_process import Activity
from indra.databases import chebi_client
import mesh_to_chembl
import requests


def process_query(drug, target,
                  save_json_name = 'chembl_output.json'):

    drug_id_db = drug.db_refs.keys()[0]

    if drug_id_db == 'CHEBI':
        chebi_id = drug.db_refs['CHEBI']
        drug_chembl_id = chebi_client.get_chembl_id(chebi_id)
    elif drug_id_db == 'MESH':
        mesh_id = drug.db_refs['MESH']
        drug_chembl_id = mesh_to_chembl.get_chembl_id(mesh_id)

    target_uid = target.db_refs['UP']

    target_chembl_id = get_target_chemblid(target_uid)

    print drug_chembl_id, target_chembl_id
          
    
    dict = chembl_client.query(drug_chembl_id, target_chembl_id)
    return process_json(dict)


def process_json(dict):
    assays = []
    for assay in dict['activities']:
        p = Activity(assay)
        p.pmid(assay)
        assays.append(p)
    return assays


def get_target_chemblid(target_uid):
    url = 'https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession=%s' % target_uid

    r = requests.get(url)
    r.raise_for_status()
    js = r.json()
    target_chemblid = js['targets'][1]['target_chembl_id']
    return target_chemblid



            
            
