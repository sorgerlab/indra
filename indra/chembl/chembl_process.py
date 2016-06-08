import requests
import json
import pandas as pd

# df = pd.read_csv('chembl_pubmed_IDs.csv',sep='\t')

class Activity(object):

    def __init__(self, assay):
        self.description = assay['assay_description']
        self.metric = assay['standard_type'] + ': ' +\
                      str(assay['standard_value']) +\
                      str(assay['standard_units'])

    def pmid(self, assay):
        chembl_doc_id = str(assay['document_chembl_id'])
        print type(chembl_doc_id)
        # pubmed_id = df[df.CHEMBL_ID == chembl_doc_id].PUBMED_ID.values[0]
        # self.pmid = 'PMID' + str(pubmed_id)
        
        url_pmid = 'https://www.ebi.ac.uk/chembl/api/data/document.json?document_chembl_id=%s' % chembl_doc_id
        r = requests.get(url_pmid)
        r.raise_for_status()
        js = r.json()
        self.pmid = 'PMID' + str(js['pubmed_id'])
        print js['pubmed_id']

        
    def __repr__(self):
        return '<%s reported in %s>' % (self.metric, self.pmid)




