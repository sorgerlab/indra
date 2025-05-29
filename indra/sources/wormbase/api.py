__all__ = ['process_from_files', 'process_from_web']


from .processor import WormBaseProcessor
from collections import namedtuple
import pandas as pd
import os

# Url for all C. elegans molecular interactions data file
alliance_mol_int_file_url = ('https://fms.alliancegenome.org/download/'
                         'INTERACTION-MOL_WB.tsv.gz')
# Url for all C. elegans genetic interactions data file
alliance_gen_int_file_url = ('https://fms.alliancegenome.org/download/'
                         'INTERACTION-GEN_WB.tsv.gz')
# Url for wormbase-to-entrez ID mapping
wormbase_entrez_mappings_file_url = ('https://ftp.ncbi.nih.gov/gene/'
                                     'DATA/GENE_INFO/Invertebrates/'
                                     'Caenorhabditis_elegans.gene_info.gz')

wormbase_int_file_url = ('https://downloads.wormbase.org/species/c_elegans/'
                         'annotation/interactions/'
                         'c_elegans.PRJNA13758.current.interactions.txt.gz')

# wormbase_all_genes_file_path = ("/Users/bradleybuchner/Desktop/Grad School/Research"
#                                 "/Aging Project/indra/indra/resources/c_elegans_all_genes.txt")
dir = os.path.dirname(os.path.abspath(__file__))
resources_dir = os.path.dirname(os.path.dirname(dir))
wormbase_all_genes_file_path = os.path.join(resources_dir, 'resources/c_elegans_all_genes.txt')


# An explanation for each column of the interaction files are here:
# https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md
alliance_int_columns = ['ids_interactor_a', 'ids_interactor_b',
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

wormbase_int_columns_raw = [
    "WBInteractionID", "Interaction_type", "Interaction_subtype", "Summary",
    "Citation", "Interactor1", "Common_name1", "Role1",
    "Interactor2", "Common_name2", "Role2",
    "Interactor3", "Common_name3", "Role3",
    "Interactor4", "Common_name4", "Role4",
    "Interactor5", "Common_name5", "Role5",
    "Interactor6", "Common_name6", "Role6",
    "Interactor7", "Common_name7", "Role7",
    "Interactor8", "Common_name8", "Role8",
    "Interactor9", "Common_name9", "Role9",
    "Interactor10", "Common_name10", "Role10",
    "Interactor11", "Common_name11", "Role11",
    "Interactor12", "Common_name12", "Role12",
    "Interactor13", "Common_name13", "Role13",
    "Interactor14", "Common_name14", "Role14",
    "Interactor15", "Common_name15", "Role15",
    "Interactor16", "Common_name16", "Role16",
    "Interactor17", "Common_name17", "Role17",
    "Interactor18", "Common_name18", "Role18",
    "Interactor19", "Common_name19", "Role19",
    "Interactor20", "Common_name20", "Role20",
    "Interactor21", "Common_name21", "Role21",
    "Interactor22", "Common_name22", "Role22",
    "Interactor23", "Common_name23", "Role23",
    "Interactor24", "Common_name24", "Role24",
    "Interactor25", "Common_name25", "Role25",
    "Interactor26", "Common_name26", "Role26",
    "Interactor27", "Common_name27", "Role27",
    "Interactor28", "Common_name28", "Role28",
    "Interactor29", "Common_name29", "Role29",
    "Interactor30", "Common_name30", "Role30",
    "Interactor31", "Common_name31", "Role31",
    "Interactor32", "Common_name32", "Role32"
    # May need to add more column names if the dataset is updated
    # with interactions having more than 32 interactors
]

wormbase_int_columns = ["WBInteractionID", "Type", "Subtype",
                                  "Effector_ID", "Effector", "Effector_Role",
                                  "Affected_ID", "Affected", "Affected_Role",
                                  "Direction"]

mapping_columns = ['tax_id', 'GeneID', 'Symbol',
                   'LocusTag', 'Synonyms', 'dbXrefs', 'chromosome',
                   'map_location', 'description', 'type_of_gene',
                   'Symbol_from_nomenclature_authority',
                   'Full_name_from_nomenclature_authority',
                   'Nomenclature_status', 'Other_designations',
                   'Modification_date', 'Feature_type']

wormbase_all_genes_columns = ['source', 'bioentity_internal_id', 'bioentity_label',
                              'synonym', 'database_xref', 'type', 'taxon',
                              'taxon_label']

_AllianceGenomeRow = namedtuple('AllianceGenomeRow', alliance_int_columns)
_WormBaseRow = namedtuple('WormBaseRow', wormbase_int_columns)


def process_from_files(alliance_gen_data_file, alliance_mol_data_file,
                       wb_all_int_data_file, wb_to_entrez_mappings_file):
    """Process WormBase interaction data from TSV files.

    Parameters
    ----------
    alliance_gen_data_file : str
        Path to the Alliance Genome dataset of C. elegans genetic interactions
        in TSV format.
    alliance_mol_data_file : str
        Path to the Alliance Genome dataset of C. elegans molecular interactions
        in TSV format.
    wb_all_int_data_file : str
        Path to the WormBase dataset of all C. elegans interactions
    wb_to_entrez_mappings_file : str
        Path to the WormBase-to-Entrez ID mapping file in TSV format.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the
        interactions data.
    """
    gen_int_iter = pd.read_csv(alliance_gen_data_file, sep='\t', comment='#',
                           dtype=str).values.tolist()
    mol_int_iter = pd.read_csv(alliance_mol_data_file, sep='\t', comment='#',
                           dtype=str).values.tolist()

    all_wb_int = pd.read_csv(wb_all_int_data_file, sep='\t', comment='#',
                             header=None, names=wormbase_int_columns_raw,
                             dtype=str, compression='gzip')
    all_wb_int_iter = _process_wormbase_interactions(all_wb_int).values.tolist()

    mappings_df = pd.read_csv(wb_to_entrez_mappings_file, sep='\t',
                              comment='#', dtype=str, names=mapping_columns)

    # return _processor_from_data(gen_int_iter, mol_int_iter, mappings_df)
    return _processor_from_data(gen_int_iter, mol_int_iter, all_wb_int_iter, mappings_df)


def process_from_web():
    """Process WormBase interaction data from the web.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the interactions data.
    """
    gen_int_iter = pd.read_csv(alliance_gen_int_file_url, sep='\t', comment='#',
                           dtype=str).values.tolist()
    mol_int_iter = pd.read_csv(alliance_mol_int_file_url, sep='\t', comment='#',
                           dtype=str).values.tolist()

    all_wb_int = pd.read_csv(wormbase_int_file_url, sep='\t', comment='#',
                             header=None, names=wormbase_int_columns_raw,
                             dtype=str, compression='gzip')
    all_wb_int_iter = _process_wormbase_interactions(all_wb_int).values.tolist()

    entrez_mappings_df = pd.read_csv(wormbase_entrez_mappings_file_url, sep='\t',
                              comment='#', dtype=str,
                              names=mapping_columns)

    all_genes_df = pd.read_csv(wormbase_all_genes_file_path, sep="\t",
                               header=None, dtype=str,
                               names=wormbase_all_genes_columns)

    # return _processor_from_data(gen_int_iter, mol_int_iter, entrez_mappings_df)
    return _processor_from_data(gen_int_iter, mol_int_iter, all_wb_int_iter,
                                entrez_mappings_df, all_genes_df)


def _process_wormbase_interactions(df):
    """Process raw interaction data file from WormBase. Expand all interactions
    into two-way interactions.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw interaction data file from WormBase.

    Returns
    -------
    processed_df : pandas.DataFrame
        Processed (pairwise) interaction data file from WormBase.
    """

    from itertools import product

    # Get relevant columns
    interactor_cols = [col for col in df.columns if col.startswith('Interactor')]
    common_name_cols = [col for col in df.columns if col.startswith('Common_name')]
    role_cols = [col for col in df.columns if col.startswith('Role')]

    # Create one row per effector-affected pair
    seen = set()
    table_rows = []
    for _, row in df.iterrows():
        participants = []
        effectors = []
        affected = []
        non_directional = []
        gene_ids = {}
        effector_roles = {}
        affected_roles = {}

        for i in range(len(interactor_cols)):
            gene_id = row.get(interactor_cols[i])
            gene_name = row.get(common_name_cols[i])
            role_raw = row.get(role_cols[i])

            if pd.isnull(gene_id) or pd.isnull(role_raw):
                continue  # skip incomplete entries

            role = str(role_raw).strip().lower()
            gene_display = gene_name if pd.notnull(gene_name) else gene_id
            gene_display = str(gene_display).strip()

            if not gene_display:
                continue  # skip blank names

            gene_ids[gene_display] = gene_id
            participants.append(gene_display)

            if role in ["effector", "target", "trans_regulator"]:
                effectors.append(gene_display)
                effector_roles[gene_display] = role
            elif role in ["affected", "bait", "trans_regulated"]:
                affected.append(gene_display)
                affected_roles[gene_display] = role
            elif role in ["non_directional"]:
                non_directional.append(gene_display)

        interaction_id = row.get("WBInteractionID", "-")
        interaction_type = row.get("Interaction_type", "-").lower()
        interaction_subtype = row.get("Interaction_subtype", "-").lower()

        if len(affected) > 0 and len(effectors) < 1:
            continue

        if len(effectors) > 0 and len(affected) < 1:
            continue

        if len(participants) < 2:
            continue

        # For directional edges
        for eff, aff in product(effectors, affected):
            key = (eff, aff, "Effector->Affected")
            if key not in seen:
                seen.add(key)
                table_rows.append({
                    "WBInteractionID": interaction_id,
                    "Type": interaction_type,
                    "Subtype": interaction_subtype,
                    "Effector_ID": gene_ids.get(eff, "-"),
                    "Effector": eff,
                    "Effector_Role": effector_roles.get(eff, "-"),
                    "Affected_ID": gene_ids.get(aff, "-"),
                    "Affected": aff,
                    "Affected_Role": affected_roles.get(aff, "-"),
                    "Direction": "Effector->Affected"
                })

        # For non-directional edges
        if not effectors and not affected and len(non_directional) >= 2:
            for i in range(len(non_directional)):
                for j in range(i + 1, len(non_directional)):
                    eff, aff = non_directional[i], non_directional[j]
                    key = tuple(sorted([eff, aff])) + ("non-directional",)
                    if key not in seen and eff != aff:
                        seen.add(key)
                        table_rows.append({
                            "WBInteractionID": interaction_id,
                            "Type": interaction_type,
                            "Subtype": interaction_subtype,
                            "Effector_ID": gene_ids.get(eff, "-"),
                            "Effector": eff,
                            "Effector_Role": "non-directional",
                            "Affected_ID": gene_ids.get(aff, "-"),
                            "Affected": aff,
                            "Affected_Role": "non-directional",
                            "Direction": "non-directional"
                        })

    # Convert to DataFrame
    processed_df = pd.DataFrame(table_rows)

    return processed_df


def _processor_from_data(gen_int_iter, mol_int_iter, all_wb_int_iter,
                         entrez_mappings_df, all_genes_df):
    """Create a WormBaseProcessor from the interaction data and ID mappings.

    Parameters
    ----------
    gen_int_iter : list
        Iterable of rows in the Alliance Genome genetic interactions data file.
    mol_int_iter : list
        Iterable of rows in the Alliance Genome molecular interactions data file.
    all_wb_int_iter : list
        Iterable of rows in the WormBase interactions data file.
    entrez_mappings_df : pd.DataFrame
        DataFrame containing associated WormBase and Entrez IDs.
    all_genes_df : pd.DataFrame
        DataFrame of all WormBase genes and their WormBase IDs,
        symbols, and synonyms.

    Returns
    -------
    indra.sources.wormbase.WormBaseProcessor
        WormBaseProcessor containing Statements extracted from the interactions data.
    """
    # Process into a list of WormBaseRow namedtuples
    alliance_rows = gen_int_iter + mol_int_iter
    alliance_data = [_AllianceGenomeRow(*[None if item == '-' else item
                           for item in row][:len(alliance_int_columns)])
            for row in alliance_rows]
    wormbase_rows = all_wb_int_iter
    wormbase_data = [_WormBaseRow(*[None if item == '-' else item
                           for item in row][:len(wormbase_int_columns)])
            for row in wormbase_rows]

    # return WormBaseProcessor(data, entrez_mappings_df)
    return WormBaseProcessor(alliance_data, wormbase_data,
                             entrez_mappings_df, all_genes_df)

