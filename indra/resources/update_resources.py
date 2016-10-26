import os
import gzip
import pandas
import rdflib
import urllib
import logging
import requests
from indra.preassembler.make_cellular_component_hierarchy import \
    get_cellular_components 

path = os.path.dirname(__file__)
logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('update_resources')
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('requests').setLevel(logging.ERROR)
logger.setLevel(logging.INFO)

def save_from_http(url, fname):
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(res.content)

def update_hgnc_entries():
    logger.info('--Updating HGNC entries-----')
    url = 'http://tinyurl.com/gnv32vh'
    fname = os.path.join(path, 'hgnc_entries.tsv')
    save_from_http(url, fname)

def update_kinases():
    logger.info('--Updating kinase list------')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=entry_name&desc=no&compress=no&query=database:(type:' + \
        'interpro%20ipr000719)%20AND%20reviewed:yes%20AND%20organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&fil=&force=no' + \
        '&format=tab&columns=id,genes(PREFERRED),organism-id,entry%20name'
    fname = os.path.join(path, 'kinases.tsv')
    save_from_http(url, fname)

def update_uniprot_entries():
    logger.info('--Updating UniProt entries--')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'fil=&force=no&format=tab&columns=id,genes(PREFERRED),organism-id,' + \
        'entry%20name'
    fname = os.path.join(path, 'uniprot_entries.tsv')
    save_from_http(url, fname)

def update_uniprot_sec_ac():
    logger.info('--Updating UniProt secondary accession--')
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
        'docs/sec_ac.txt'
    logger.info('Downloading %s' % url)
    fname = os.path.join(path, 'uniprot_sec_ac.txt')
    urllib.urlretrieve(url, fname)

def update_uniprot_subcell_loc():
    logger.info('--Updating UniProt subcellular location--')
    url = 'http://www.uniprot.org/locations/?' + \
        '%20sort=&desc=&compress=no&query=&force=no&format=tab&columns=id'
    fname = os.path.join(path, 'uniprot_subcell_loc.tsv')
    save_from_http(url, fname)

def update_chebi_entries():
    logger.info('--Updating ChEBI entries----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/reference.tsv.gz'
    fname = os.path.join(path, 'reference.tsv.gz')
    urllib.urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
    # Save PubChem mapping
    fname = os.path.join(path, 'chebi_to_pubchem.tsv')
    logger.info('Saving into %s' % fname)
    df_pubchem = df[df['REFERENCE_DB_NAME']=='PubChem']
    df_pubchem.sort_values('COMPOUND_ID', ascending=True, inplace=True)
    df_pubchem.to_csv(fname, sep='\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                      header=['CHEBI', 'PUBCHEM'], index=False)
    # Save ChEMBL mapping
    fname = os.path.join(path, 'chebi_to_chembl.tsv')
    logger.info('Saving into %s' % fname)
    df_chembl = df[df['REFERENCE_DB_NAME']=='ChEMBL']
    df_chembl.sort_values('COMPOUND_ID', ascending=True, inplace=True)
    df_chembl.to_csv(fname, sep='\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                      header=['CHEBI', 'CHEMBL'], index=False)

def update_cellular_components():
    logger.info('--Updating GO cellular components----')
    url = 'http://purl.obolibrary.org/obo/go.owl'
    fname = os.path.join(path, 'go.owl')
    save_from_http(url, fname)
    g = rdflib.Graph()
    g.parse(fname)
    component_map, component_part_map = get_cellular_components(g)
    fname = os.path.join(path, 'cellular_components.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write('id\tname\n')

        for comp_id, comp_name in sorted(component_map.items(),
                                          key=lambda x: x[0]):
            fh.write('%s\t%s\n' % (comp_id, comp_name))

if __name__ == '__main__':
    update_hgnc_entries()
    update_kinases()
    update_uniprot_entries()
    update_uniprot_sec_ac()
    update_uniprot_subcell_loc()
    update_chebi_entries()
    update_cellular_components()
