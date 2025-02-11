from .processor import WormBaseProcessor
import requests
import gzip
from io import BytesIO
from collections import namedtuple
import pandas as pd


wormbase_mol_file_url = ('https://fms.alliancegenome.org/download/'
                     'INTERACTION-MOL_WB.tsv.gz') # Url for all C. elegans molecular interactions data file
wormbase_gen_file_url = ('https://fms.alliancegenome.org/download/'
                     'INTERACTION-GEN_WB.tsv.gz') # Url for all C. elegans genetic interactions data file
wormbase_entrez_mappings_file_url = ('https://ftp.ncbi.nih.gov/gene/'
                                     'DATA/GENE_INFO/Invertebrates/'
                                     'Caenorhabditis_elegans.gene_info.gz') # Url for wormbase-to-entrez ID mapping

# An explanation for each column of the interaction files are here:
# https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md
columns = ['ids_interactor_a', 'ids_interactor_b',
           'alt_ids_interactor_a', 'alt_ids_interactor_b',
           'aliases_interactor_a', 'aliases_interactor_b',
           'interaction_detection_methods', 'publication_first_authors',
           'publication_identifiers', 'taxid_interactor_a',
           'taxid_interactor_b', 'interaction_types',
           'source_databases', 'interaction_identifiers',
           'confidence_values', 'expansion_methods',
           'biological_roles_interactor_a',
           'biological_roles_interactor_b',
           'experimental_roles_interactor_a',
           'experimental_roles_interactor_b',
           'types_interactor_a', 'types_interactor_b',
           'xrefs_interactor_a', 'xrefs_interactor_b',
           'interaction_xrefs', 'annotations_interactor_a',
           'annotations_interactor_b', 'interaction_annotations',
           'host_organisms', 'interaction_parameters',
           'creation_date', 'update_date', 'checksums_interactor_a',
           'checksums_interactor_b', 'interaction_checksums',
           'negative', 'features_interactor_a', 'features_interactor_b',
           'stoichiometries_interactor_a', 'stoichiometries_interactor_b',
           'identification_method_participant_a',
           'identification_method_participant_b']

mapping_columns = ['tax_id', 'GeneID', 'Symbol',
                   'LocusTag', 'Synonyms', 'dbXrefs', 'chromosome',
                   'map_location', 'description', 'type_of_gene', 'Symbol_from_nomenclature_authority',
                   'Full_name_from_nomenclature_authority',	'Nomenclature_status', 'Other_designations',
                   'Modification_date', 'Feature_type']

_WormBaseRow = namedtuple('WormBaseRow', columns)

def process_from_file(wormbase_gen_data_file, wormbase_mol_data_file,
                      wb_to_entrez_mappings_file):
    """Process WormBase interaction data from TSV files.

    Parameters
    ----------
    wormbase_gen_data_file : str
        Path to the WormBase genetic interactions data file in TSV format.
    wormbase_mol_data_file : str
        Path to the WormBase molecular interactions data file in TSV format.
    wb_to_entrez_mappings_file : str
        Path to the WormBase-to-Entrez ID mapping file in TSV format.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the interactions data.
    """
    gen_iter = pd.read_csv(wormbase_gen_data_file, sep='\t', comment='#', na_values="-", dtype=str).values.tolist()
    mol_iter = pd.read_csv(wormbase_mol_data_file, sep='\t', comment='#', na_values="-", dtype=str).values.tolist()
    mappings_df = pd.read_csv(wb_to_entrez_mappings_file, sep='\t', comment='#', dtype=str, names=mapping_columns)
    return _processor_from_data(gen_iter, mol_iter, mappings_df)


def process_from_web(wormbase_gen_data_file=None, wormbase_mol_data_file=None,
                     wb_to_entrez_mappings_file=None):
    """Process WormBase interaction data from the web.

    Parameters
    ----------
    wormbase_gen_data_file : Optional[str]
        If specified, the genetic interactions data will be written to this file.
    wormbase_mol_data_file : Optional[str]
        If specified, the molecular interactions data will be written to this file.
    wb_to_entrez_mappings_file : Optional[str]
        If specified, the WormBase-to-Entrez ID mapping file will be written to this file.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the interactions data.
    """

    res_gen = requests.get(wormbase_gen_file_url)
    if res_gen.status_code == 200:
        gzip_bytes_gen = BytesIO(res_gen.content)
        with gzip.open(gzip_bytes_gen, 'rt') as file:
            gen_iter = pd.read_csv(file, sep='\t', comment='#', na_values="-", dtype=str).values.tolist()
    else:
        raise Exception('Unable to download WormBase genetic interactions data: status code %s'
                        % res_gen.status_code)

    res_mol = requests.get(wormbase_mol_file_url)
    if res_mol.status_code == 200:
        gzip_bytes_mol = BytesIO(res_mol.content)
        with gzip.open(gzip_bytes_mol, 'rt') as file:
            mol_iter = pd.read_csv(file, sep='\t', comment='#', na_values="-", dtype=str).values.tolist()
    else:
        raise Exception('Unable to download WormBase molecular interactions data: status code %s'
                        % res_mol.status_code)

    res_map = requests.get(wormbase_entrez_mappings_file_url)
    if res_map.status_code == 200:
        gzip_bytes_map = BytesIO(res_map.content)
        with gzip.open(gzip_bytes_map, 'rt') as file:
            mappings_df = pd.read_csv(file, sep='\t', comment='#', na_values="-", dtype=str, names=mapping_columns)
    else:
        raise Exception('Unable to download WormBase-to-Entrez ID mapping data: status code %s'
                        % res_mol.status_code)

    return _processor_from_data(gen_iter, mol_iter, mappings_df)


def _processor_from_data(gen_iter, mol_iter, mappings_df):
    """Create a WormBaseProcessor from the interaction data and ID mappings.

    Parameters
    ----------
    gen_iter : list
        Iterable of rows in the genetic interactions data file.
    mol_iter : list
        Iterable of rows in the molecular interactions data file.
    mappings_df : pd.DataFrame
        DataFrame containing associated WormBase and Entrez IDs.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the interactions data.
    """
    # Process into a list of WormBaseRow namedtuples
    all_rows = gen_iter + mol_iter
    data = [_WormBaseRow(*[None if item == '-' else item
                            for item in row][:len(columns)]) for row in all_rows]
    return WormBaseProcessor(data, mappings_df)
