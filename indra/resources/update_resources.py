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

def load_from_http(url):
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    return res.content

def save_from_http(url, fname):
    content = load_from_http(url)
    if content is None:
        return
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(content)

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
    reviewed_entries = load_from_http(url)
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=&fil=organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&force=no&' + \
        'format=tab&columns=id,genes(PREFERRED),organism-id,entry%20name'
    unreviewed_human_entries = load_from_http(url)
    if not((reviewed_entries is not None) and
            (unreviewed_human_entries is not None)):
            return
    lines = reviewed_entries.strip().split('\n')
    lines += unreviewed_human_entries.strip().split('\n')[1:]
    full_table = '\n'.join(lines)
    fname = os.path.join(path, 'uniprot_entries.tsv')
    with open(fname, 'wb') as fh:
        fh.write(full_table)

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

def update_bel_chebi_map():
    logger.info('--Updating BEL ChEBI map----')
    url = 'http://resource.belframework.org/belframework/latest-release/' + \
        'equivalence/'
    url1 = url + 'chebi-ids.beleq'
    url2 = url + 'chebi.beleq'
    res1 = load_from_http(url1)
    res2 = load_from_http(url2)
    if not ((res1 is not None) and (res2 is not None)):
        return
    id_lines = [lin.strip() for lin in res1.split('\n')]
    started = False
    id_map = {}
    for id_line in id_lines:
        if started:
            if id_line:
                # Instead of splitting on |, split using UUID fixed length
                chebi_id = id_line[:-37]
                uuid = id_line[-36:]
                id_map[uuid] = chebi_id
        if id_line == '[Values]':
            started = True
    name_lines = [lin.strip() for lin in res2.split('\n')]
    started = False
    name_map = {}
    for name_line in name_lines:
        if started:
            if name_line:
                # Instead of splitting on |, split using UUID fixed length
                chebi_name = name_line[:-37]
                uuid = name_line[-36:]
                name_map[uuid] = chebi_name
        if name_line == '[Values]':
            started = True
    fname = os.path.join(path, 'bel_chebi_map.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        for uuid, chebi_name in sorted(name_map.items(), key=lambda x: x[1]):
            chebi_id = id_map.get(uuid)
            if chebi_id is not None:
                fh.write('%s\tCHEBI:%s\n' % (chebi_name, chebi_id))

if __name__ == '__main__':
    #update_hgnc_entries()
    #update_kinases()
    update_uniprot_entries()
    #update_uniprot_sec_ac()
    #update_uniprot_subcell_loc()
    #update_chebi_entries()
    #update_cellular_components()
    #update_bel_chebi_map()
