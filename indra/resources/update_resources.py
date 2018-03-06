from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import gzip
import pandas
import rdflib
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve
import logging
import requests
from indra.util import read_unicode_csv, write_unicode_csv
from indra.preassembler.make_cellular_component_hierarchy import \
    get_cellular_components
from indra.preassembler.make_cellular_component_hierarchy import \
    main as make_ccomp_hierarchy
from indra.preassembler.make_entity_hierarchy import \
    main as make_ent_hierarchy
from indra.preassembler.make_activity_hierarchy import \
    main as make_act_hierarchy
from indra.preassembler.make_modification_hierarchy import \
    main as make_mod_hierarchy

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
    url = 'http://tinyurl.com/y83dx5s6'
    fname = os.path.join(path, 'hgnc_entries.tsv')
    save_from_http(url, fname)

def update_kinases():
    logger.info('--Updating kinase list------')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=entry_name&desc=no&compress=no&query=database:(type:' + \
        'interpro%20ipr011009)%20AND%20reviewed:yes%20AND%20organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&fil=&force=no' + \
        '&format=tab&columns=id,genes(PREFERRED),organism-id,entry%20name'
    fname = os.path.join(path, 'kinases.tsv')
    save_from_http(url, fname)

def update_uniprot_entries():
    logger.info('--Updating UniProt entries--')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'format=tab&columns=id,genes(PREFERRED),' + \
        'entry%20name,database(RGD),database(MGI)'
    reviewed_entries = load_from_http(url)
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:no&fil=organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&' + \
        'format=tab&columns=id,genes(PREFERRED),entry%20name,' + \
        'database(RGD),database(MGI)'
    unreviewed_human_entries = load_from_http(url)
    if not((reviewed_entries is not None) and
            (unreviewed_human_entries is not None)):
            return
    lines = reviewed_entries.strip('\n').split('\n')
    lines += unreviewed_human_entries.strip('\n').split('\n')[1:]
    # At this point, we need to clean up the gene names.
    logging.info('Processing UniProt entries list.')
    for i, line in enumerate(lines):
        if i == 0:
            continue
        terms = line.split('\t')
        # If there are multiple gene names, take the first one
        gene_names = terms[1].split(';')
        terms[1] = gene_names[0]
        # Join the line again after the change
        lines[i] = '\t'.join(terms)
    # Join all lines into a single string
    full_table = '\n'.join(lines)
    fname = os.path.join(path, 'uniprot_entries.tsv')
    logging.info('Saving into %s.' % fname)
    with open(fname, 'wb') as fh:
        fh.write(full_table)

def update_uniprot_sec_ac():
    logger.info('--Updating UniProt secondary accession--')
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
        'docs/sec_ac.txt'
    logger.info('Downloading %s' % url)
    fname = os.path.join(path, 'uniprot_sec_ac.txt')
    urlretrieve(url, fname)

def update_uniprot_subcell_loc():
    # TODO: This file could be stored as a tsv instead after some processing
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
    urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
    # Save PubChem mapping
    fname = os.path.join(path, 'chebi_to_pubchem.tsv')
    logger.info('Saving into %s' % fname)
    df_pubchem = df[df['REFERENCE_DB_NAME']=='PubChem']
    df_pubchem.sort_values(['COMPOUND_ID', 'REFERENCE_ID'], ascending=True,
                           inplace=True)
    df_pubchem.to_csv(fname, sep=b'\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                      header=['CHEBI', 'PUBCHEM'], index=False)
    # Save ChEMBL mapping
    fname = os.path.join(path, 'chebi_to_chembl.tsv')
    logger.info('Saving into %s' % fname)
    df_chembl = df[df['REFERENCE_DB_NAME']=='ChEMBL']
    df_chembl.sort_values(['COMPOUND_ID', 'REFERENCE_ID'], ascending=True,
                          inplace=True)
    df_chembl.to_csv(fname, sep=b'\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                      header=['CHEBI', 'CHEMBL'], index=False)


def update_cas_to_chebi():
    logger.info('--Updating CAS to ChEBI entries----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/database_accession.tsv'
    fname = os.path.join(path, 'database_accession.tsv')
    urlretrieve(url, fname)
    with open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
    fname = os.path.join(path, 'cas_to_chebi.tsv')
    logger.info('Saving into %s' % fname)
    df_cas = df[df['TYPE'] == 'CAS Registry Number']
    df_cas.sort_values(['ACCESSION_NUMBER', 'COMPOUND_ID'], ascending=True,
                       inplace=True)
    # Here we need to map to primary ChEBI IDs
    with open('chebi_to_primary.tsv', 'rb') as fh:
        df_prim = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
        mapping = {s: p for s, p in zip(df_prim['Secondary'].tolist(),
                                        df_prim['Primary'].tolist())}
    df_cas.COMPOUND_ID.replace(mapping, inplace=True)
    df_cas.drop_duplicates(subset=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                           inplace=True)
    df_cas.to_csv(fname, sep=b'\t',
                  columns=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                  header=['CAS', 'CHEBI'], index=False)


def update_chebi_primary_map():
    logger.info('--Updating ChEBI primary map entries----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/compounds.tsv.gz'
    fname = os.path.join(path, 'compounds.tsv.gz')
    urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
    fname = os.path.join(path, 'chebi_to_primary.tsv')
    logger.info('Saving into %s' % fname)
    df = df[df['PARENT_ID'] != 'null']
    df.replace('CHEBI:([0-9]+)', r'\1', inplace=True, regex=True)
    df.sort_values(['CHEBI_ACCESSION', 'PARENT_ID'], ascending=True,
                   inplace=True)
    df.drop_duplicates(subset=['CHEBI_ACCESSION', 'PARENT_ID'], inplace=True)
    df.to_csv(fname, sep=b'\t',
              columns=['CHEBI_ACCESSION', 'PARENT_ID'], 
              header=['Secondary', 'Primary'], index=False)


def update_cellular_components():
    logger.info('--Updating GO cellular components----')
    url = 'http://purl.obolibrary.org/obo/go.owl'
    fname = os.path.join(path, '../../data/go.owl')
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
    id_lines = []
    name_lines = []
    for release in ('1.0', 'latest-release'):
        url = 'http://resource.belframework.org/belframework/%s/' % release + \
              'equivalence/'
        url1 = url + 'chebi-ids.beleq'
        if release == '1.0':
            url2 = url + 'chebi-names.beleq'
        else:
            url2 = url + 'chebi.beleq'
        res1 = load_from_http(url1)
        res2 = load_from_http(url2)
        id_lines1 = [lin.strip() for lin in res1.split('\n') if lin]
        start = id_lines1.index('[Values]')
        id_lines1 = id_lines1[start+1:]
        id_lines += id_lines1
        name_lines1 = [lin.strip() for lin in res2.split('\n') if lin]
        start = name_lines1.index('[Values]')
        name_lines1 = name_lines1[start + 1:]
        name_lines += name_lines1
    id_map = {}
    for id_line in id_lines:
        if id_line:
            # Instead of splitting on |, split using UUID fixed length
            chebi_id = id_line[:-37]
            uuid = id_line[-36:]
            id_map[uuid] = chebi_id

    name_map = {}
    for name_line in name_lines:
        if name_line:
            # Instead of splitting on |, split using UUID fixed length
            chebi_name = name_line[:-37]
            uuid = name_line[-36:]
            name_map[uuid] = chebi_name

    name_to_id = {}
    for uuid, chebi_name in name_map.items():
        chebi_id = id_map.get(uuid)
        if chebi_id is not None:
            if chebi_name in name_to_id:
                old_id = int(name_to_id[chebi_name])
                new_id = int(chebi_id)
                if new_id <= old_id:
                    continue
            name_to_id[chebi_name] = chebi_id

    fname = os.path.join(path, 'bel_chebi_map.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        for chebi_name, chebi_id in sorted(name_to_id.items(),
                                           key=lambda x: x[0]):
            fh.write('%s\tCHEBI:%s\n' % (chebi_name, chebi_id))

def update_entity_hierarchy():
    logger.info('--Updating entity hierarchy----')
    fname = os.path.join(path, 'famplex/relations.csv')
    make_ent_hierarchy(fname)

def update_modification_hierarchy():
    logger.info('--Updating modification hierarchy----')
    make_mod_hierarchy()

def update_activity_hierarchy():
    logger.info('--Updating activity hierarchy----')
    make_act_hierarchy()

def update_cellular_component_hierarchy():
    logger.info('--Updating cellular component hierarchy----')
    make_ccomp_hierarchy()

def update_famplex_map():
    logger.info('--Updating FamPlex map----')
    # Currently this is a trivial "copy" of the FamPlex equivalences.csv
    # file. Later, name spaces may need to be adapted and other format changes
    # may be needed.
    fname_in = os.path.join(path, 'famplex/equivalences.csv')
    fname_out = os.path.join(path, 'famplex_map.tsv')
    rows = read_unicode_csv(fname_in)
    write_unicode_csv(fname_out, rows, delimiter='\t')

def update_ncit_map():
    logger.info('--Updating NCIT map----')
    url_hgnc = 'https://ncit.nci.nih.gov/ncitbrowser/ajax?action=' + \
               'export_mapping&dictionary=NCIt_to_HGNC_Mapping&version=1.0'

    url_go = 'https://ncit.nci.nih.gov/ncitbrowser/ajax?action=' + \
             'export_mapping&dictionary=GO_to_NCIt_Mapping&version=1.1'

    url_chebi = 'https://ncit.nci.nih.gov/ncitbrowser/ajax?action=' + \
                'export_mapping&dictionary=NCIt_to_ChEBI_Mapping&version=1.0'

    def get_ncit_df(url):
        df = pandas.read_csv(url)
        df = df[df['Association Name'] == 'mapsTo']
        df.sort_values(['Source Code', 'Target Code'], ascending=True,
                       inplace=True)
        df = df[['Source Code', 'Target Code', 'Source Coding Scheme',
                 'Target Coding Scheme']]
        return df

    df_hgnc = get_ncit_df(url_hgnc)
    df_hgnc.replace('HGNC:(\d*)\s*', '\\1', inplace=True, regex=True)
    df_go = get_ncit_df(url_go)
    df_go.rename(columns={'Source Code': 'Target Code',
                       'Target Code': 'Source Code',
                       'Source Coding Scheme': 'Target Coding Scheme',
                       'Target Coding Scheme': 'Source Coding Scheme'},
              inplace=True)
    df_chebi = get_ncit_df(url_chebi)
    df_chebi.replace('ChEBI', 'CHEBI', inplace=True)

    # Add the old HGNC mappings
    df_hgnc_old = pandas.read_csv('ncit_hgnc_map_old.tsv', sep='\t',
                                  index_col=None, dtype=str)
    df_hgnc = df_hgnc.append(df_hgnc_old)
    df_hgnc.sort_values(['Source Code', 'Target Code'], ascending=True,
                        inplace=True)

    # Add UniProt mappings
    df_uniprot = pandas.read_csv('Feb2017NCIt-SwissProt.txt', sep='\t',
                                 index_col=None)
    up_entries = {'Source Code': [], 'Target Coding Scheme': [],
                  'Target Code': []}
    for entry in df_uniprot.iterrows():
        up_entries['Source Code'].append(entry[1]['code'].strip())
        up_entries['Target Coding Scheme'].append('UP')
        up_entries['Target Code'].append(entry[1]['Swiss_Prot'].strip())
    df_uniprot = pandas.DataFrame.from_dict(up_entries)
    df_uniprot.sort_values(['Source Code', 'Target Code'], ascending=True,
                           inplace=True)

    df_all = pandas.concat([df_chebi, df_go, df_hgnc, df_uniprot])

    fname = os.path.join(path, 'ncit_map.tsv')
    df_all.to_csv(fname, sep=b'\t', columns=['Source Code',
                                        'Target Coding Scheme',
                                        'Target Code'],
              header=['NCIT ID', 'Target NS', 'Target ID'], index=False)

def update_chebi_names():
    logger.info('--Updating ChEBI names----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/names_3star.tsv.gz'
    fname = os.path.join(path, 'names_3star.tsv.gz')
    urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
    fname = os.path.join(path, 'chebi_names.tsv')
    df = df[df['TYPE'] == 'NAME']
    df.sort_values(by='COMPOUND_ID', inplace=True)
    logger.info('Saving into %s' % fname)
    df.to_csv(fname, sep=b'\t', header=True, index=False,
              columns=['COMPOUND_ID', 'NAME'])


def update_famplex():
    """Update all the CSV files that form the FamPlex resource."""
    famplex_url_pattern = \
        'https://raw.githubusercontent.com/sorgerlab/famplex/master/%s.csv'
    csv_names = ['entities', 'equivalences', 'gene_prefixes',
                 'grounding_map', 'relations']
    for csv_name in csv_names:
        url = famplex_url_pattern % csv_name
        save_from_http(url, 'famplex/%s.csv' % csv_name)


if __name__ == '__main__':
    update_famplex()
    update_hgnc_entries()
    update_kinases()
    update_uniprot_entries()
    update_uniprot_sec_ac()
    update_uniprot_subcell_loc()
    update_chebi_entries()
    update_chebi_names()
    update_cas_to_chebi()
    update_cellular_components()
    update_bel_chebi_map()
    update_entity_hierarchy()
    update_modification_hierarchy()
    update_activity_hierarchy()
    update_cellular_component_hierarchy()
    update_famplex_map()
    update_ncit_map()
    update_chebi_primary_map()
