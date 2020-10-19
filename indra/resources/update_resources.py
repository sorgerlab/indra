import os
import json
import gzip
import pandas
import logging
import requests
from zipfile import ZipFile
from collections import defaultdict
from urllib.request import urlretrieve
from xml.etree import ElementTree as ET
from indra.util import read_unicode_csv, write_unicode_csv
from indra.databases.obo_client import OboClient
from indra.databases import chebi_client, pubchem_client
from indra.databases.lincs_client import load_lincs_csv
from . import load_resource_json, get_resource_path

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
            'md_mgd_id', 'md_rgd_id', 'gd_prev_sym', 'gd_pub_ensembl_id']
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


def update_chebi_references():
    # The reference table contains all the automated mappings from ChEBI
    # IDs to IDs in other databases, except CAS, which only has manually
    # curated mappings available in the database_accession table
    # (see implementation in update_cas_to_chebi).
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
    # If the second column of the row starts with CID:, strip out the CID
    # prefix Otherwise, include the row unchanged
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
            # Include other rows unchanged
            new_rows.append(original_row)
    write_unicode_csv(fname, new_rows, '\t')

    # In another round of cleanup, we try dealing with duplicate mappings in a
    # principled way such that many-to-one mappings are allowed but one-to-many
    # mappings are eliminated
    original_rows = read_unicode_csv(fname, '\t')
    chebi_pubchem = defaultdict(list)
    pubchem_chebi = defaultdict(list)
    for chebi_id, pc_id in original_rows:
        chebi_pubchem[chebi_id].append(pc_id)
        pubchem_chebi[pc_id].append(chebi_id)
    # Looking for InChIKey matches for duplicates in the ChEBI -> PubChem
    # direction
    logger.info('Getting InChiKey matches for duplicates')
    ik_matches = set()
    for chebi_id, pc_ids in chebi_pubchem.items():
        if len(pc_ids) > 1:
            ck = chebi_client.get_inchi_key(chebi_id)
            for pc_id in pc_ids:
                pk = pubchem_client.get_inchi_key(pc_id)
                if ck == pk:
                    ik_matches.add((chebi_id, pc_id))
    # Looking for InChIKey matches for duplicates in the PubChem -> ChEBI
    # direction
    for pc_id, chebi_ids in pubchem_chebi.items():
        if len(chebi_ids) > 1:
            pk = pubchem_client.get_inchi_key(pc_id)
            for chebi_id in chebi_ids:
                ck = chebi_client.get_inchi_key(chebi_id)
                if ck == pk:
                    ik_matches.add((chebi_id, pc_id))
    rows = read_unicode_csv(fname, '\t')
    header = next(rows)
    header.append('IK_MATCH')
    new_rows = [header]
    for chebi_id, pc_id in rows:
        if (chebi_id, pc_id) in ik_matches:
            new_rows.append([chebi_id, pc_id, 'Y'])
        else:
            new_rows.append([chebi_id, pc_id, ''])
    write_unicode_csv(fname, new_rows, '\t')

    # Save ChEMBL mapping
    fname = os.path.join(path, 'chebi_to_chembl.tsv')
    logger.info('Saving into %s' % fname)
    df_chembl = df[df['REFERENCE_DB_NAME'] == 'ChEMBL']
    # Get additional mappings for compounds in tas
    df_chembl_tas = pandas.read_csv(os.path.join(path, 'chembl_tas.csv'),
                                    sep=',')[['chebi_id', 'chembl_id']]
    df_chembl_tas = df_chembl_tas[~df_chembl_tas.chebi_id.isna()]
    df_chembl_tas['chebi_id'] = df_chembl_tas.chebi_id.\
        apply(lambda x: str(int(x)))
    df_chembl_tas.columns = ['COMPOUND_ID', 'REFERENCE_ID']
    df_chembl = pandas.concat([df_chembl, df_chembl_tas]).drop_duplicates()
    df_chembl.sort_values(['COMPOUND_ID', 'REFERENCE_ID'], ascending=True,
                          inplace=True)
    df_chembl.to_csv(fname, sep='\t', columns=['COMPOUND_ID', 'REFERENCE_ID'],
                     header=['CHEBI', 'CHEMBL'], index=False)


def update_hmdb_chebi_map():
    logger.info('--Updating HMDB to ChEBI entries----')
    ns = {'hmdb': 'http://www.hmdb.ca'}
    url = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'
    fname = os.path.join(path, 'hmdb_metabolites.zip')
    logger.info('Downloading %s' % url)
    #urlretrieve(url, fname)
    mappings = []
    with ZipFile(fname) as input_zip:
        with input_zip.open('hmdb_metabolites.xml') as fh:
            for event, elem in ET.iterparse(fh, events=('start', 'end')):
                #print(elem.tag)
                if event == 'start' and \
                        elem.tag == '{%s}metabolite' % ns['hmdb']:
                    hmdb_id = None
                    chebi_id = None
                # Important: we only look at accession if there's no HMDB
                # ID yet, otherwise we pick up secondary accession tags
                elif event == 'start' and \
                        elem.tag == '{%s}accession' % ns['hmdb'] and \
                        not hmdb_id:
                    hmdb_id = elem.text
                elif event == 'start' and \
                        elem.tag == '{%s}chebi_id' % ns['hmdb']:
                    chebi_id = elem.text
                elif event == 'end' and \
                        elem.tag == '{%s}metabolite' % ns['hmdb']:
                    if hmdb_id and chebi_id:
                        print(hmdb_id, chebi_id)
                        mappings.append([hmdb_id, chebi_id])
                elem.clear()
    fname = os.path.join(path, 'hmdb_to_chebi.tsv')
    mappings = [['HMDB_ID', 'CHEBI_ID']] + sorted(mappings, key=lambda x: x[0])
    write_unicode_csv(fname, mappings, delimiter='\t')


def update_chebi_accessions():
    # The database_accession table contains manually curated mappings
    # between ChEBI and other databases. It only contains very few mappings
    # to e.g., PubChem, therefore the main resource for those mappings
    # is the reference table (see implementation in update_chebi_entries).
    # The only useful mappings we extract here from the database_accession
    # table are the ones to CAS which are not available in the reference
    # table.
    logger.info('--Updating CAS to ChEBI entries----')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/' + \
        'Flat_file_tab_delimited/database_accession.tsv'
    fname = os.path.join(path, 'database_accession.tsv')
    urlretrieve(url, fname)
    with open(fname, 'rb') as fh:
        logger.info('Loading %s' % fname)
        df = pandas.read_csv(fh, sep='\t', index_col=None,
                             dtype=str, na_filter=False)
    fname = os.path.join(path, 'cas_to_chebi.tsv')
    logger.info('Saving into %s' % fname)
    df_cas = df[df['TYPE'] == 'CAS Registry Number']
    df_cas.sort_values(['ACCESSION_NUMBER', 'COMPOUND_ID'], ascending=True,
                       inplace=True)
    # Here we need to map to primary ChEBI IDs
    from indra.databases.chebi_client import get_primary_id
    df_cas.COMPOUND_ID.replace(get_primary_id, inplace=True)
    df_cas.drop_duplicates(subset=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                           inplace=True)
    df_cas.to_csv(fname, sep='\t',
                  columns=['ACCESSION_NUMBER', 'COMPOUND_ID'],
                  header=['CAS', 'CHEBI'], index=False)


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
            fh.write(('%s\tCHEBI:%s\n' %
                      (chebi_name, chebi_id)).encode('utf-8'))


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

    url_swissprot = 'https://evs.nci.nih.gov/ftp1/' \
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
    allele_fname = os.path.join(path, 'ncit_allele_map.tsv')
    df_hgnc_old = pandas.read_csv(allele_fname, sep='\t',
                                  index_col=None, dtype=str)
    df_hgnc = df_hgnc.append(df_hgnc_old)
    df_hgnc.sort_values(['Source Code', 'Target Code'], ascending=True,
                        inplace=True)

    # Add UniProt mappings
    ncit_swissprot_file = 'NCIt-SwissProt.txt'
    # These entries represent non-Uniprot (nuccode?) IDs and
    # we exclude them for validity
    exclude_list = {'M74558', 'U37690', 'U65002'}
    save_from_http(url_swissprot, ncit_swissprot_file)
    df_uniprot = pandas.read_csv(ncit_swissprot_file, sep='\t',
                                 index_col=None)
    up_entries = {'Source Code': [], 'Target Coding Scheme': [],
                  'Target Code': []}
    for entry in df_uniprot.iterrows():
        up_id = entry[1]['SwissProt ID'].strip()
        if up_id in exclude_list:
            continue
        up_entries['Source Code'].append(entry[1]['NCIt Code'].strip())
        up_entries['Target Coding Scheme'].append('UP')
        up_entries['Target Code'].append(up_id)
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
        save_from_http(url, os.path.join(path, 'famplex',
                                         '%s.csv' % csv_name))


def update_grounding_map():
    famplex_gmap = os.path.join(path, 'famplex', 'grounding_map.csv')
    covid_gmap = os.path.join(path, 'grounding', 'covid_grounding.csv')
    famplex_rows = list(read_unicode_csv(famplex_gmap))
    row_len = len(famplex_rows[0])
    covid_rows = list(read_unicode_csv(covid_gmap))
    covid_rows = [r + [''] * (row_len - len(r))
                  for r in covid_rows]
    all_rows = famplex_rows + covid_rows
    grounding_map = os.path.join(path, 'grounding', 'grounding_map.csv')
    write_unicode_csv(grounding_map, all_rows)


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
    """Update Mesh ID to name and tree number mappings."""
    url = ('ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/'
           'xmlmesh/desc2021.gz')
    desc_path = os.path.join(path, 'mesh_desc2021.gz')
    if not os.path.exists(desc_path):
        logging.info('Download MeSH descriptors from %s', url)
        urlretrieve(url, desc_path)
        logging.info('Done downloading MeSH descriptors')
    # Process the XML and find descriptor records
    with gzip.open(desc_path) as desc_file:
        logging.info('Parsing MeSH descriptors')
        et = ET.parse(desc_file)
    rows = []
    for record in et.iterfind('DescriptorRecord'):
        # We first get the ID and the name
        uid = record.find('DescriptorUI').text
        name = record.find('DescriptorName/String').text
        term_name_str = _get_term_name_str(record, name)
        tree_numbers = record.findall('TreeNumberList/TreeNumber')
        tree_numbers_str = '|'.join([t.text for t in tree_numbers])
        rows.append((uid, name, term_name_str, tree_numbers_str))

    fname = os.path.join(path, 'mesh_id_label_mappings.tsv')
    write_unicode_csv(fname, rows, delimiter='\t')


def update_mesh_supplementary_names():
    """Update MeSH ID to name mappings for supplementary terms."""
    supp_url = ('ftp://nlmpubs.nlm.nih.gov/online/mesh/MESH_FILES/'
                'xmlmesh/supp2021.gz')
    supp_path = os.path.join(path, 'mesh_supp2021.gz')
    if not os.path.exists(supp_path):
        logging.info('Download MeSH supplement from %s', supp_url)
        urlretrieve(supp_url, supp_path)
        logging.info('Done downloading MeSH supplement')
    with gzip.open(supp_path) as supp_file:
        logging.info('Parsing MeSH supplement')
        supp_et = ET.parse(supp_file)
    supp_rows = []
    for record in supp_et.iterfind('SupplementalRecord'):
        uid = record.find('SupplementalRecordUI').text
        name = record.find('SupplementalRecordName/String').text
        term_name_str = _get_term_name_str(record, name)
        supp_rows.append((uid, name, term_name_str))

    fname = os.path.join(path, 'mesh_supp_id_label_mappings.tsv')
    write_unicode_csv(fname, supp_rows, delimiter='\t')


def _get_term_name_str(record, name):
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
    return term_name_str


def update_mirbase():
    url = 'ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz'

    mirbase_gz_path = os.path.join(path, 'mirbase.gz')
    if not os.path.exists(mirbase_gz_path):
        urlretrieve(url, mirbase_gz_path)

    mirbase_json_path = os.path.join(path, 'mirbase.tsv')
    with gzip.open(mirbase_gz_path, 'rt') as in_file, \
              open(mirbase_json_path, 'w') as out_file:
        print('mirbase_id', 'mirbase_name', 'database', 'identifier', 'label',
              sep='\t', file=out_file)
        for line in _process_mirbase_file(in_file):
            print(*line, sep='\t', file=out_file)

    # This is a big intermediate file, so don't keep it
    os.remove(mirbase_gz_path)


def _process_mirbase_file(lines):
    """Process the lines of the miRBase definitions file."""
    groups = []

    for line in lines:
        if line.startswith('ID'):
            groups.append([])
        groups[-1].append(line)

    for group in groups:
        mirbase_name = group[0][5:23].strip()
        mirbase_id = group[2][3:-2].strip()
        for element in group:
            if not element.startswith('DR'):
                continue
            db, identifier, name = [e.strip() for e in \
                                    element[len('DR'):].lstrip().split(';')]
            yield mirbase_id, mirbase_name, db, identifier, name.rstrip('.')


def update_doid():
    """Update disease ontology."""
    url = 'http://purl.obolibrary.org/obo/doid.obo'
    OboClient.update_resource(path, url, 'doid', remove_prefix=False)


def update_efo():
    """Update experimental factor ontology."""
    url = 'https://www.ebi.ac.uk/efo/efo.obo'
    OboClient.update_resource(path, url, 'efo', remove_prefix=True)


def update_hpo():
    """Update human phenotype ontology."""
    url = 'http://purl.obolibrary.org/obo/hp.obo'
    OboClient.update_resource(path, url, 'hp', remove_prefix=False)


def update_go():
    """Update gene ontology."""
    url = 'http://purl.obolibrary.org/obo/go.obo'
    OboClient.update_resource(path, url, 'go', remove_prefix=False)


def update_chebi_obo():
    """Update disease ontology."""
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo'
    OboClient.update_resource(path, url, 'chebi', remove_prefix=False)


def update_drugbank_mappings():
    """Update mappings from DrugBank to CHEBI/CHEMBL"""
    # Note that for this to work, PyOBO (https://github.com/pyobo/pyobo) has
    # to be installed and the DrugBank download
    # (https://www.drugbank.ca/releases/latest) put into ~/.obo/drugbank/
    # Note that the DrugBank download requires signing up for an account and
    # waiting for approval.
    import pyobo
    drugbank_chembl = pyobo.get_filtered_xrefs('drugbank', 'chembl.compound')
    drugbank_chebi = pyobo.get_filtered_xrefs('drugbank', 'chebi')
    chebi_drugbank = pyobo.get_filtered_xrefs('chebi', 'drugbank')
    drugbank_names = pyobo.get_id_name_mapping('drugbank')
    rows = []
    for drugbank_id, chembl_id in drugbank_chembl.items():
        rows.append([drugbank_id, 'CHEMBL', chembl_id, 'drugbank'])
    for drugbank_id, chebi_id in drugbank_chebi.items():
        rows.append([drugbank_id, 'CHEBI', chebi_id, 'drugbank'])
    for chebi_id, drugbank_id in chebi_drugbank.items():
        rows.append([drugbank_id, 'CHEBI', chebi_id, 'chebi'])
    for drugbank_id, name in drugbank_names.items():
        rows.append([drugbank_id, 'NAME', name, 'drugbank'])
    fname = os.path.join(path, 'drugbank_mappings.tsv')
    header = ['DRUGBANK_ID', 'NAMESPACE', 'ID', 'SOURCE']
    rows = [header] + sorted(rows)
    write_unicode_csv(fname, rows, delimiter='\t')


def update_identifiers_registry():
    """Update prefixes and patterns for identifiers namespaces."""
    url = \
        'https://registry.api.identifiers.org/resolutionApi/getResolverDataset'
    res = requests.get(url)
    regj = res.json()
    patterns = {entry['prefix']:
                    {'pattern': entry['pattern'],
                     'namespace_embedded': entry['namespaceEmbeddedInLui']}
                for entry in sorted(regj['payload']['namespaces'],
                                    key=lambda x: x['prefix'])}
    with open(os.path.join(path, 'identifiers_patterns.json'), 'w') as fh:
        json.dump(patterns, fh, indent=1)


def update_biomappings():
    """Update mappings from the BioMappings project."""
    from indra.databases import mesh_client
    from indra.databases.identifiers import get_ns_id_from_identifiers
    biomappings_branch = 'import_gilda_mappings'

    # Predictions
    url = (f'https://raw.githubusercontent.com/biomappings/biomappings/'
           f'{biomappings_branch}/predictions.tsv')
    df_pred = pandas.read_csv(url, delimiter='\t', dtype=str)
    # Mappings
    url = (f'https://raw.githubusercontent.com/biomappings/biomappings/'
           f'{biomappings_branch}/mappings.tsv')
    df_mappings = pandas.read_csv(url, delimiter='\t', dtype=str)

    df = df_mappings.append(df_pred)

    ### Apply any necessary filters here
    # We take only exact matches
    df = df[df['relation'] == 'skos:exactMatch']
    # We only take real xrefs, not refs within a given ontology
    df = df[df['source prefix'] != df['target prefix']]

    # We now construct a mapping dict of these mappings
    biomappings = {}
    for _, row in df.iterrows():
        source_ns, source_id = \
            get_ns_id_from_identifiers(row['source prefix'],
                                       row['source identifier'])
        target_ns, target_id = \
            get_ns_id_from_identifiers(row['target prefix'],
                                       row['target identifier'])
        biomappings[(source_ns, source_id, row['source name'])] = \
            (target_ns, target_id, row['target name'])
        biomappings[(target_ns, target_id, row['target name'])] = \
            (source_ns, source_id, row['target name'])

    mesh_mappings = {k: v for k, v in biomappings.items() if k[0] == 'MESH'}
    non_mesh_mappings = {k: v for k, v in biomappings.items()
                         if k[0] != 'MESH' and v[0] != 'MESH'}
    rows = [list(k + v) for k, v in
            sorted(non_mesh_mappings.items(), key=lambda x: x[0][1])]
    write_unicode_csv(get_resource_path('biomappings.tsv'), rows,
                      delimiter='\t')

    # We next look at mappings to MeSH from EFO/HP/DOID
    for ns in ['efo', 'hp', 'doid']:
        for entry in load_resource_json('%s.json' % ns):
            db, db_id, name = ns.upper(), entry['id'], entry['name']
            if (db, db_id) in biomappings:
                continue
            # We first need to decide if we prioritize another name space
            xref_dict = {xr['namespace']: xr['id']
                         for xr in entry.get('xrefs', [])}
            if 'MESH' in xref_dict or 'MSH' in xref_dict:
                mesh_id = xref_dict.get('MESH') or xref_dict.get('MSH')
                if not mesh_id.startswith('D'):
                    continue
                mesh_name = mesh_client.get_mesh_name(mesh_id)
                if not mesh_name:
                    continue
                mesh_mappings[('MESH', mesh_id, name)] = \
                    (db, db_id, entry['name'])
    rows = [list(k + v) for k, v in
            sorted(mesh_mappings.items(), key=lambda x: x[0][1])]
    write_unicode_csv(get_resource_path('mesh_mappings.tsv'), rows,
                      delimiter='\t')


def main():
    update_famplex()
    update_famplex_map()
    update_grounding_map()
    update_hgnc_entries()
    update_kinases()
    update_uniprot_subcell_loc()
    update_chebi_obo()
    update_chebi_references()
    update_chebi_accessions()
    update_hmdb_chebi_map()
    update_bel_chebi_map()
    update_ncit_map()
    update_lincs_small_molecules()
    update_lincs_proteins()
    update_mesh_names()
    update_mesh_supplementary_names()
    update_mirbase()
    update_doid()
    update_efo()
    update_hpo()
    update_drugbank_mappings()
    update_identifiers_registry()


if __name__ == '__main__':
    main()
