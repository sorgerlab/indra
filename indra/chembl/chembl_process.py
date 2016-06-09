import requests
import json


class Activity(object):

    def __init__(self, assay):
        self.description = assay['assay_description']
        self.metric = assay['standard_type'] + ': ' +\
                      str(assay['standard_value']) +\
                      str(assay['standard_units'])

    def pmid(self, assay):
        chembl_doc_id = str(assay['document_chembl_id'])      
        url_pmid = 'https://www.ebi.ac.uk/chembl/api/data/document.json?document_chembl_id=%s' % chembl_doc_id
        r = requests.get(url_pmid)
        r.raise_for_status()
        js = r.json()
        self.pmid = 'PMID' + str(js['documents'][0]['pubmed_id'])
       
        
    def __repr__(self):
        return '<%s reported in %s>' % (self.metric, self.pmid)




