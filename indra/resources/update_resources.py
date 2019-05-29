from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import gzip
import pickle
import pandas
import rdflib
from indra.util import read_unicode_csv, write_unicode_csv

try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve
import logging
import requests
from indra.util import read_unicode_csv, write_unicode_csv
from indra.databases import go_client
from indra.databases import uniprot_client
from indra.databases.lincs_client import load_lincs_csv
from indra.preassembler import make_cellular_component_hierarchy as mcch
from indra.preassembler.make_entity_hierarchy import \
    main as make_ent_hierarchy
from indra.preassembler.make_activity_hierarchy import \
    main as make_act_hierarchy
from indra.preassembler.make_modification_hierarchy import \
    main as make_mod_hierarchy
from indra.preassembler.hierarchy_manager import get_bio_hierarchies

path = os.path.dirname(__file__)
logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('indra.resources.update_resources')
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('requests').setLevel(logging.ERROR)
logger.setLevel(logging.INFO)

# Set a global variable indicating whether we've downloaded the latest GO
# during this update cycle so that we don't do it more than once
go_updated = False


def load_latest_go():
    global go_updated
    go_fname = go_client.go_owl_path
    if not go_updated:
        url = 'http://purl.obolibrary.org/obo/go.owl'
        print("Downloading latest GO from %s" % url)
        save_from_http(url, go_fname)
        go_updated = True
    return go_client.load_go_graph(go_fname)


def load_from_http(url):
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download "%s"' % url)
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

    # Select relevant columns and parameters
    cols = ['gd_hgnc_id', 'gd_app_sym', 'gd_app_name', 'gd_status',
            'gd_aliases', 'md_eg_id', 'md_prot_id',
            'md_mgd_id', 'md_rgd_id', 'gd_prev_sym']
    statuses = ['Approved', 'Entry%20Withdrawn']
    params = {
            'hgnc_dbtag': 'on',
            'order_by': 'gd_app_sym_sort',
            'format': 'text',
            'submit': 'submit'
            }

    # Construct a download URL from the above parameters
    url = 'https://www.genenames.org/cgi-bin/download/custom?'
    url += '&'.join(['col=%s' % c for c in cols]) + '&'
    url += '&'.join(['status=%s' % s for s in statuses]) + '&'
    url += '&'.join(['%s=%s' % (k, v) for k, v in params.items()])

    # Save the download into a file
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

    from indra.databases import hgnc_client, uniprot_client
    add_kinases = ['PGK1', 'PKM', 'TAF1', 'NME1', 'BCKDK', 'PDK1', 'PDK2',
                   'PDK3', 'PDK4', 'BCR', 'FAM20C', 'BAZ1B', 'PIKFYVE']
    df = pandas.read_csv(fname, sep='\t')
    for kinase in add_kinases:
        hgnc_id = hgnc_client.get_hgnc_id(kinase)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        up_mnemonic = uniprot_client.get_mnemonic(up_id)
        df = df.append({'Entry': up_id, 'Gene names  (primary )': kinase,
                        'Organism ID': '9606', 'Entry name': up_mnemonic},
                       ignore_index=True)
    df.to_csv(fname, sep='\t', index=False)


def update_uniprot_subcell_loc():
    logger.info('--Updating UniProt subcellular location--')
    url = 'https://www.uniprot.org/docs/subcell.txt'
    res = requests.get(url)
    res.raise_for_status()
    header, entry_block = res.text.split('_' * 75)
    entries = entry_block.split('//')
    mappings = []
    for entry in entries:
        slid = None
        goid = None
        lines = entry.split('\n')
        for line in lines:
            if line.startswith('AC'):
                slid = line[5:].strip()
            if line.startswith('GO'):
                goid = line[5:].split(';')[0]
        if slid and goid:
            mappings.append((slid, goid))
    fname = os.path.join(path, 'uniprot_subcell_loc.tsv')
    write_unicode_csv(fname, mappings, delimiter='\t')


def update_chebi_entries():
    logger.info('--Updating ChEBI entries----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/reference.tsv.gz'
    fname = os.path.join(path, 'reference.tsv.gz')
    urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.read_csv(fh, sep='\t', index_col=None,
                             parse_dates=True, encoding='latin-1')
    # Save PubChem mapping
    fname = os.path.join(path, 'chebi_to_pubchem.tsv')
    logger.info('Saving into %s' % fname)
    df_pubchem = df[df['REFERENCE_DB_NAME']=='PubChem']
    df_pubchem.sort_values(['COMPOUND_ID', 'REFERENCE_ID'], ascending=True,
                           inplace=True)
    df_pubchem.to_csv(fname, sep='\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                      header=['CHEBI', 'PUBCHEM'], index=False)

    # Process PubChem mapping to eliminate SID rows and strip CID: prefix
    # If the second column of the row starts with SID:, ignore the row
    # If the second column of the row starts with CID:, strip out the CID prefix
    # Otherwise, include the row unchanged
    original_rows = read_unicode_csv(fname, '\t')
    new_rows = []
    for original_row in original_rows:
        if original_row[1].startswith('CID:'):
            new_row = original_row
            new_row[1] = new_row[1][5:] # Strip out CID:
            new_rows.append(new_row)
        elif original_row[1].startswith('SID:'):
            # Skip SID rows
            continue
        else:
            # Include other rows unchanges
            new_rows.append(original_row)
    write_unicode_csv(fname, new_rows, '\t')

    # Save ChEMBL mapping
    fname = os.path.join(path, 'chebi_to_chembl.tsv')
    logger.info('Saving into %s' % fname)
    df_chembl = df[df['REFERENCE_DB_NAME']=='ChEMBL']
    df_chembl.sort_values(['COMPOUND_ID', 'REFERENCE_ID'], ascending=True,
                          inplace=True)
    df_chembl.to_csv(fname, sep='\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
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
    with open(os.path.join(path, 'chebi_to_primary.tsv'), 'rb') as fh:
        df_prim = pandas.DataFrame.from_csv(fh, sep='\t', index_col=None)
        mapping = {s: p for s, p in zip(df_prim['Secondary'].tolist(),
                                        df_prim['Primary'].tolist())}
    df_cas.COMPOUND_ID.replace(mapping, inplace=True)
    df_cas.drop_duplicates(subset=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                           inplace=True)
    df_cas.to_csv(fname, sep='\t',
                  columns=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                  header=['CAS', 'CHEBI'], index=False)


def update_chebi_primary_map_and_names():
    logger.info('--Updating ChEBI primary map entries and names----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/compounds.tsv.gz'
    fname = os.path.join(path, 'compounds.tsv.gz')
    urlretrieve(url, fname)
    with gzip.open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df_orig = pandas.read_csv(fh, sep='\t', index_col=None,
                                  parse_dates=True, dtype='str')
    # This df is still shared
    df_orig.replace('CHEBI:([0-9]+)', r'\1', inplace=True, regex=True)

    # First we construct mappings to primary accession IDs
    # This df is specific to parents
    df = df_orig[df_orig['PARENT_ID'].notna()]
    df.sort_values(['CHEBI_ACCESSION', 'PARENT_ID'], ascending=True,
                   inplace=True)
    df.drop_duplicates(subset=['CHEBI_ACCESSION', 'PARENT_ID'], inplace=True)
    fname = os.path.join(path, 'chebi_to_primary.tsv')
    logger.info('Saving into %s' % fname)
    df.to_csv(fname, sep='\t',
              columns=['CHEBI_ACCESSION', 'PARENT_ID'], 
              header=['Secondary', 'Primary'], index=False)

    # Second we get the ID to name mappings
    df = df_orig[df_orig['NAME'].notna()]
    df.sort_values(by=['CHEBI_ACCESSION', 'ID'], inplace=True)

    fname = os.path.join(path, 'chebi_names.tsv')
    logger.info('Saving into %s' % fname)
    df.to_csv(fname, sep='\t', header=True, index=False,
              columns=['CHEBI_ACCESSION', 'NAME'])


def update_cellular_component_hierarchy():
    logger.info('--Updating GO cellular components----')
    g = load_latest_go()
    component_map, component_part_map = go_client.get_cellular_components(g)
    # Save the cellular component ID->name mappings
    fname = os.path.join(path, 'cellular_components.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write('id\tname\n'.encode('utf-8'))
        for comp_id, comp_name in sorted(component_map.items(),
                                          key=lambda x: x[0]):
            fh.write(('%s\t%s\n' % (comp_id, comp_name)).encode('utf-8'))
    # Create the cellular component hierarchy
    gg = mcch.make_component_hierarchy(component_map, component_part_map)
    mcch.save_hierarchy(gg, mcch.rdf_file)


def update_go_id_mappings():
    g = load_latest_go()
    go_client.update_id_mappings(g)


def update_bel_chebi_map():
    logger.info('--Updating BEL ChEBI map----')
    id_lines = []
    name_lines = []

    url = 'https://raw.githubusercontent.com/OpenBEL/' + \
            'openbel-framework-resources/latest/equivalence/'
    url1 = url + 'chebi-ids.beleq'
    url2 = url + 'chebi.beleq'
    res1 = load_from_http(url1).decode('utf-8')
    res2 = load_from_http(url2).decode('utf-8')
    id_lines1 = [lin.strip() for lin in res1.split('\n') if lin]
    start = id_lines1.index('[Values]')
    id_lines1 = id_lines1[start+1:]
    id_lines += id_lines1
    name_lines1 = [lin.strip() for lin in res2.split('\n') if lin]
    start = name_lines1.index('[Values]')
    name_lines1 = name_lines1[start + 1:]
    name_lines += name_lines1

    # Here we need to get the ChEBI to primary map to make sure we map
    # to primary IDs
    with open(os.path.join(path, 'chebi_to_primary.tsv'), 'r') as fh:
        chebi_to_primary = {k: v for k, v in
                            [l.strip().split('\t') for
                             l in fh.readlines()][1:]}

    id_map = {}
    for id_line in id_lines:
        if id_line:
            # Instead of splitting on |, split using UUID fixed length
            chebi_id = id_line[:-37]
            uuid = id_line[-36:]
            # Map secondary IDs to primary IDs before adding to the map
            if chebi_id in chebi_to_primary:
                chebi_id = chebi_to_primary[chebi_id]
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
            fh.write(('%s\tCHEBI:%s\n' % (chebi_name, chebi_id)).encode('utf-8'))


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

    url_swissprot = 'https://ncit.nci.nih.gov/ncitbrowser/ajax?action=' \
                    'export_mapping&uri=https://evs.nci.nih.gov/ftp1/' \
                    'NCI_Thesaurus/Mappings/NCIt-SwissProt_Mapping.txt'

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
    df_hgnc_old = pandas.read_csv('ncit_allele_map.tsv', sep='\t',
                                  index_col=None, dtype=str)
    df_hgnc = df_hgnc.append(df_hgnc_old)
    df_hgnc.sort_values(['Source Code', 'Target Code'], ascending=True,
                        inplace=True)

    # Add UniProt mappings
    ncit_swissprot_file = 'NCIt-SwissProt.txt'
    save_from_http(url_swissprot, ncit_swissprot_file)
    df_uniprot = pandas.read_csv(ncit_swissprot_file, sep='\t',
                                 index_col=None)
    up_entries = {'Source Code': [], 'Target Coding Scheme': [],
                  'Target Code': []}
    for entry in df_uniprot.iterrows():
        up_entries['Source Code'].append(entry[1]['NCIt Code'].strip())
        up_entries['Target Coding Scheme'].append('UP')
        up_entries['Target Code'].append(entry[1]['SwissProt ID'].strip())
    df_uniprot = pandas.DataFrame.from_dict(up_entries)
    df_uniprot.sort_values(['Source Code', 'Target Code'], ascending=True,
                           inplace=True)

    df_all = pandas.concat([df_chebi, df_go, df_hgnc, df_uniprot])

    fname = os.path.join(path, 'ncit_map.tsv')
    df_all.to_csv(fname, sep='\t', columns=['Source Code',
                                            'Target Coding Scheme',
                                            'Target Code'],
                  header=['NCIT ID', 'Target NS', 'Target ID'], index=False)


def update_famplex():
    """Update all the CSV files that form the FamPlex resource."""
    famplex_url_pattern = \
        'https://raw.githubusercontent.com/sorgerlab/famplex/master/%s.csv'
    csv_names = ['entities', 'equivalences', 'gene_prefixes',
                 'grounding_map', 'relations']
    for csv_name in csv_names:
        url = famplex_url_pattern % csv_name
        save_from_http(url, os.path.join(path,'famplex/%s.csv' % csv_name))


def update_hierarchy_pickle():
    fname = os.path.join(path, 'bio_hierarchies.pkl')
    hierarchies = get_bio_hierarchies(from_pickle=False)
    with open(fname, 'wb') as fh:
        pickle.dump(hierarchies, fh, protocol=4)


def update_lincs_small_molecules():
    """Load the csv of LINCS small molecule metadata into a dict.

    Produces a dict keyed by HMS LINCS small molecule ids, with the metadata
    contained in a dict of row values keyed by the column headers extracted
    from the csv.
    """
    url = 'http://lincs.hms.harvard.edu/db/sm/'
    sm_data = load_lincs_csv(url)
    sm_dict = {d['HMS LINCS ID']: d.copy() for d in sm_data}
    assert len(sm_dict) == len(sm_data), "We lost data."
    fname = os.path.join(path, 'lincs_small_molecules.json')
    with open(fname, 'w') as fh:
        json.dump(sm_dict, fh, indent=1)


def update_lincs_proteins():
    """Load the csv of LINCS protein metadata into a dict.

    Produces a dict keyed by HMS LINCS protein ids, with the metadata
    contained in a dict of row values keyed by the column headers extracted
    from the csv.
    """
    url = 'http://lincs.hms.harvard.edu/db/proteins/'
    prot_data = load_lincs_csv(url)
    prot_dict = {d['HMS LINCS ID']: d.copy() for d in prot_data}
    assert len(prot_dict) == len(prot_data), "We lost data."
    fname = os.path.join(path, 'lincs_proteins.json')
    with open(fname, 'w') as fh:
        json.dump(prot_dict, fh, indent=1)


def update_mesh_names():
    url = 'ftp://nlmpubs.nlm.nih.gov/online/mesh/2018/xmlmesh/desc2018.xml'
    urlretrieve(url, 'desc2018.xml')
    # Process the XML and find descriptor records
    et = ET.parse('desc2018.xml')
    records = et.findall('DescriptorRecord')
    rows = []
    for record in records:
        # We first get the ID and the name
        uid = record.find('DescriptorUI').text
        name = record.find('DescriptorName/String').text
        # We then need to look for additional terms related to the
        # preferred concept to get additional names
        concepts = record.findall('ConceptList/Concept')
        all_term_names = []
        for concept in concepts:
            # We only look at the preferred concept here
            if concept.attrib['PreferredConceptYN'] == 'Y':
                terms = concept.findall('TermList/Term')
                for term in terms:
                    term_name = term.find('String').text
                    if term_name != name:
                        all_term_names.append(term_name)
        # Append a list of term names separated by pipes to the table
        term_name_str = '|'.join(all_term_names)
        rows.append((uid, name, term_name_str))
    fname = os.path.join(path, 'mesh_id_label_mappings.tsv')
    write_unicode_csv(fname, rows, delimiter='\t')


if __name__ == '__main__':
    update_go_id_mappings()
    update_cellular_component_hierarchy()
    update_famplex()
    update_famplex_map()
    update_hgnc_entries()
    update_kinases()
    update_uniprot_subcell_loc()
    update_chebi_entries()
    update_chebi_primary_map_and_names()
    update_cas_to_chebi()
    update_bel_chebi_map()
    update_entity_hierarchy()
    update_modification_hierarchy()
    update_activity_hierarchy()
    update_hierarchy_pickle()
    update_ncit_map()
    update_lincs_small_molecules()
    update_lincs_proteins()
    update_mesh_names()
