import requests
import json


def query(drug_chemblid, target_chemblid):
    
    url = 'https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=%s&target_chembl_id=%s' % (drug_chemblid, target_chemblid)

    r = requests.get(url)
    r.raise_for_status()
    js = r.json()
    return js

