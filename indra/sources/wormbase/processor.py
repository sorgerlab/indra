
import re
import csv
import os
import tqdm
import logging
import requests
import pandas as pd
import numpy as np
from io import BytesIO
import gzip
from collections import namedtuple
from indra.statements import Agent, Evidence, Activation, Association, Inhibition, Phosphorylation, Demethylation, \
    Methylation
from indra.ontology.standardize import standardize_name_db_refs


logger = logging.getLogger(__name__)

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

_WormBaseRow = namedtuple('WormBaseRow', columns)

"""Miscellaneous info for WormBase interaction data (genetic and molecular):

    Unique source databases: 
    ['wormbase' 'biogrid' 'MINT' 'IntAct' 'UniProt' 'DIP']

    Unique agent ID types: 
    ['wormbase' 'entrez gene/locuslink' 'uniprotkb' 'intact'] 

    Unique interaction ID types: 
    ['wormbase' 'biogrid' 'intact' 'mint' 'imex' 'dip' 'wwpdb' 'emdb']

"""

class WormBaseProcessor(object):
    """Extracts INDRA statements from WormBase interaction data.

        Parameters
        ----------
        wormbase_gen_file : str
            The file from WormBase containing all genetic interactions in
            C. elegans. If not provided, the WormBase data is downloaded
            from the WormBase website (technically alliancegenome.org).
        wormbase_mol_file : str
            The file from WormBase containing all molecular interactions in
            C. elegans. If not provided, the WormBase data is downloaded
            from the WormBase website (technically alliancegenome.org).

        Attributes
        ----------
        statements : list[indra.statements.Statements]
            Extracted INDRA statements.
        """

    def __init__(self, wormbase_gen_file=None, wormbase_mol_file=None):
        self.statements = []
        self.wormbase_gen_file = wormbase_gen_file
        self.wormbase_mol_file = wormbase_mol_file

        # Load the wormbase-to-entrez ID mapping file as a pd.DataFrame
        wb_to_entrez_mappings = self._download_wormbase_data(wormbase_entrez_mappings_file_url,
                                                             as_df=True)

        # If a path to the genetic interactions file is provided, process it.
        if self.wormbase_gen_file:
            rows_gen = self._read_wormbase_data(self.wormbase_gen_file)
        # If no file is provided, download from web
        else:
            logger.info('No genetic interactions file specified, downloading from WormBase '
                        'at %s' % wormbase_gen_file_url)
            rows_gen = self._download_wormbase_data(wormbase_gen_file_url, as_df=False)

        # If a path to the molecular interactions file is provided, process it.
        if self.wormbase_mol_file:
            rows_mol = self._read_wormbase_data(self.wormbase_mol_file)
        # If no file is provided, download from web
        else:
            logger.info('No molecular interactions file specified, downloading from WormBase '
                        'at %s' % wormbase_mol_file_url)
            rows_mol = self._download_wormbase_data(wormbase_mol_file_url, as_df=False)

        rows = rows_gen + rows_mol

        # Convert 'dbXrefs' column in wb_to_entrez_mappings to a dictionary of Source:ID pairs.
        wb_to_entrez_mappings['dbXrefs'] = wb_to_entrez_mappings['dbXrefs'].apply(self._id_conversion)

        # Create new column 'wormbase_id' that holds each gene's WormBase identifier
        wb_to_entrez_mappings['wormbase_id'] = wb_to_entrez_mappings['dbXrefs'].apply(
            lambda x: x.get('WormBase')[0] if isinstance(x, dict) and 'WormBase' in x
            else x.get('WB')[0] if isinstance(x, dict) and 'WB' in x
            else None
        )

        # Convert mappings to dictionaries for quick lookups
        wb_to_entrez_dict = wb_to_entrez_mappings.set_index('wormbase_id')['GeneID'].to_dict()
        entrez_to_wb_dict = wb_to_entrez_mappings.set_index('GeneID')['wormbase_id'].to_dict()
        entrez_to_symbol_dict = wb_to_entrez_mappings.set_index('GeneID')['Symbol'].to_dict()
        symbol_to_annotation_dict = wb_to_entrez_mappings.set_index('GeneID').to_dict(orient='index')

        # Process the rows into Statements
        for idx, row in enumerate(tqdm.tqdm(rows, desc='Processing WormBase rows')):
            try:
                # There are some extra columns that we don't need to take and
                # thereby save space in annotations
                filt_row = [None if item == '-' else item
                            for item in row][:len(columns)]
                wb_row = _WormBaseRow(*filt_row)

                # Get the name of agent A
                name_agent_a = None
                alias_info_agent_a = self._alias_conversion(wb_row.aliases_interactor_a) or {}
                alt_ids_agent_a = self._id_conversion(wb_row.alt_ids_interactor_a) or {}
                if not alias_info_agent_a: # If agent alias is empty, look for a valid name in alternate IDs
                    if not alt_ids_agent_a:
                        logger.warning(f"Agent alias and alternate ID dicts for interactor A are empty: {wb_row}")
                    else: # If the alternate ids dict is not empty, look for names in the order below, with 'entrez
                        # gene/locuslink' and lowercase preferred.
                        all_lowercase_names = []
                        all_uppercase_names = []
                        for key in ['entrez gene/locuslink', 'uniprot/swiss-prot', 'biogrid']:
                            if alt_ids_agent_a.get(key):
                                lowercase_names = [s for s in (alt_ids_agent_a.get(key) or []) if s.islower()]
                                uppercase_names = [s for s in (alt_ids_agent_a.get(key) or []) if not s.islower()]
                                if lowercase_names:
                                    all_lowercase_names.extend(lowercase_names)
                                if uppercase_names:
                                    all_uppercase_names.extend(uppercase_names)
                        if all_lowercase_names:
                            name_agent_a = all_lowercase_names[0]
                        elif all_uppercase_names:
                            name_agent_a = all_uppercase_names[0]
                        else:
                            # If no names were found above, use whatever first value is in the alt. ids dict
                            # as a fallback
                            name_agent_a = next(iter(alt_ids_agent_a.values()), [None])[0]
                            # logger.warning(f"No valid lowercase or uppercase name found for interactor A. Using "
                            #                f"fallback: {name_agent_a} ... ... {wb_row}")
                else: # If the alias dict is not empty, look for names in the order below, with 'public_name' and
                    # lowercase preferred.
                    all_lowercase_names = []
                    all_uppercase_names = []
                    for key in ['public_name', 'gene name', 'display_short', 'gene name synonym']:
                        if alias_info_agent_a.get(key):
                            lowercase_names = [s for s in (alias_info_agent_a.get(key) or []) if s.islower()]
                            uppercase_names = [s for s in (alias_info_agent_a.get(key) or []) if not s.islower()]
                            if lowercase_names:
                                all_lowercase_names.extend(lowercase_names)
                            if uppercase_names:
                                all_uppercase_names.extend(uppercase_names)
                    if all_lowercase_names:
                        name_agent_a = all_lowercase_names[0]
                    elif all_uppercase_names:
                        name_agent_a = all_uppercase_names[0]
                    else:
                        # If no names were found above, use whatever first value is in the alias dict
                        # as a fallback
                        name_agent_a = next(iter(alias_info_agent_a.values()), [None])[0]
                        # logger.warning(f"No valid lowercase or uppercase name found for interactor A. Using "
                        #                f"fallback: {name_agent_a} ... ... {wb_row}")

                # Get the name of agent B
                name_agent_b = None
                alias_info_agent_b = self._alias_conversion(wb_row.aliases_interactor_b) or {}
                alt_ids_agent_b = self._id_conversion(wb_row.alt_ids_interactor_b) or {}
                if not alias_info_agent_b:  # If agent alias is empty, look for a valid name in alternate IDs
                    if not alt_ids_agent_b:
                        logger.warning(
                            f"Agent alias and alternate ID dicts for interactor B are empty: {wb_row}")
                    else:  # If the alternate ids dict is not empty, look for names in the order below, with 'entrez
                        # gene/locuslink' and lowercase preferred.
                        all_lowercase_names = []
                        all_uppercase_names = []
                        for key in ['entrez gene/locuslink', 'uniprot/swiss-prot', 'biogrid']:
                            if alt_ids_agent_b.get(key):
                                lowercase_names = [s for s in (alt_ids_agent_b.get(key) or []) if s.islower()]
                                uppercase_names = [s for s in (alt_ids_agent_b.get(key) or []) if
                                                   not s.islower()]
                                if lowercase_names:
                                    all_lowercase_names.extend(lowercase_names)
                                if uppercase_names:
                                    all_uppercase_names.extend(uppercase_names)
                        if all_lowercase_names:
                            name_agent_b = all_lowercase_names[0]
                        elif all_uppercase_names:
                            name_agent_b = all_uppercase_names[0]
                        else:
                            # If no names were found above, use whatever first value is in the alt. ids dict
                            # as a fallback
                            name_agent_b = next(iter(alt_ids_agent_b.values()), [None])[0]
                            # logger.warning(
                            #     f"No valid lowercase or uppercase name found for interactor B. Using "
                            #     f"fallback: {name_agent_b} ... ... {wb_row}")
                else:  # If the alias dict is not empty, look for names in the order below, with 'public_name' and
                    # lowercase preferred.
                    all_lowercase_names = []
                    all_uppercase_names = []
                    for key in ['public_name', 'gene name', 'display_short', 'gene name synonym']:
                        if alias_info_agent_b.get(key):
                            lowercase_names = [s for s in (alias_info_agent_b.get(key) or []) if s.islower()]
                            uppercase_names = [s for s in (alias_info_agent_b.get(key) or []) if
                                               not s.islower()]
                            if lowercase_names:
                                all_lowercase_names.extend(lowercase_names)
                            if uppercase_names:
                                all_uppercase_names.extend(uppercase_names)
                    if all_lowercase_names:
                        name_agent_b = all_lowercase_names[0]
                    elif all_uppercase_names:
                        name_agent_b = all_uppercase_names[0]
                    else:
                        # If no names were found above, use whatever first value is in the alias dict
                        # as a fallback
                        name_agent_b = next(iter(alias_info_agent_b.values()), [None])[0]
                        # logger.warning(f"No valid lowercase or uppercase name found for interactor B. Using "
                        #                f"fallback: {name_agent_b} ... ... {wb_row}")


                # Get db_refs using wb_row.ids_interactor_(a/b)
                wormbase_id_agent_a = None
                entrez_id_agent_a = None
                up_id_agent_a = None
                intact_id_agent_a = None
                db_id_info_agent_a = self._id_conversion(wb_row.ids_interactor_a) or {}
                alt_db_id_info_agent_a = self._id_conversion(wb_row.alt_ids_interactor_a) or {}

                if not db_id_info_agent_a:
                    logger.warning(f"No db_refs found for interactor A: {wb_row}")
                else:
                    if db_id_info_agent_a.get('wormbase'):
                        wormbase_id_agent_a = db_id_info_agent_a.get('wormbase')[0]
                    # Some WB ids are stored as an alternate id under 'ensemblgenomes'
                    elif alt_db_id_info_agent_a.get('ensemblgenomes') and 'WBGene' in \
                        alt_db_id_info_agent_a.get('ensemblgenomes'):
                            wormbase_id_agent_a = alt_db_id_info_agent_a.get('ensemblgenomes')[0]
                    if db_id_info_agent_a.get('entrez gene/locuslink'):
                        entrez_id_agent_a = db_id_info_agent_a.get('entrez gene/locuslink')[0]
                    elif wormbase_id_agent_a: # If an entrez ID isn't found but a WB ID is, use mappings file to get
                        entrez_id_agent_a = wb_to_entrez_dict.get(wormbase_id_agent_a) or None

                    if not wormbase_id_agent_a and entrez_id_agent_a: # If WB ID isn't found but an entrez ID is,
                        wormbase_id_agent_a = entrez_to_wb_dict.get(entrez_id_agent_a) or None

                    if db_id_info_agent_a.get('uniprotkb'):
                        up_id_agent_a = db_id_info_agent_a.get('uniprotkb')[0]
                    if db_id_info_agent_a.get('intact'):
                        intact_id_agent_a = db_id_info_agent_a.get('intact')[0]

                wormbase_id_agent_b = None
                entrez_id_agent_b = None
                up_id_agent_b = None
                intact_id_agent_b = None
                db_id_info_agent_b = self._id_conversion(wb_row.ids_interactor_b) or {}
                alt_db_id_info_agent_b = self._id_conversion(wb_row.alt_ids_interactor_b) or {}
                if not db_id_info_agent_b:
                    logger.warning(f"No db_refs found for interactor B: {wb_row}")
                else:
                    if db_id_info_agent_b.get('wormbase'):
                        wormbase_id_agent_b = db_id_info_agent_b.get('wormbase')[0]
                    # Some WB ids are stored as an alternate id under 'ensemblgenomes'
                    elif alt_db_id_info_agent_b.get('ensemblgenomes') and 'WBGene' in \
                         alt_db_id_info_agent_b.get('ensemblgenomes'):
                        wormbase_id_agent_b = alt_db_id_info_agent_b.get('ensemblgenomes')[0]
                    if db_id_info_agent_b.get('entrez gene/locuslink'):
                        entrez_id_agent_b = db_id_info_agent_b.get('entrez gene/locuslink')[0]
                    elif wormbase_id_agent_b: # If an entrez ID isn't found but a WB ID is, use mappings file to get
                        entrez_id_agent_b = wb_to_entrez_dict.get(wormbase_id_agent_b) or None

                    if not wormbase_id_agent_b and entrez_id_agent_b: # If WB ID isn't found but an entrez ID is,
                        wormbase_id_agent_b = entrez_to_wb_dict.get(entrez_id_agent_b) or None

                    if db_id_info_agent_b.get('uniprotkb'):
                        up_id_agent_b = db_id_info_agent_b.get('uniprotkb')[0]
                    if db_id_info_agent_b.get('intact'):
                        intact_id_agent_b = db_id_info_agent_b.get('intact')[0]

                # If agent name doesn't match the corresponding name in the wormbase-to-entrez ID mapping file, replace
                # it with the name in that file.
                if entrez_id_agent_a:
                    entrez_name_agent_a = entrez_to_symbol_dict.get(entrez_id_agent_a) or None
                    if entrez_name_agent_a and name_agent_a and name_agent_a != entrez_name_agent_a:
                        logger.warning(f"Replacing name for interactor A with Entrez symbol: {name_agent_a} "
                                       f"--> {entrez_name_agent_a}")
                        name_agent_a = entrez_name_agent_a

                if entrez_id_agent_b:
                    entrez_name_agent_b = entrez_to_symbol_dict.get(entrez_id_agent_b) or None
                    if entrez_name_agent_b and name_agent_b and name_agent_b != entrez_name_agent_b:
                        logger.warning(f"Replacing name for interactor B with Entrez symbol: {name_agent_b} "
                                       f"--> {entrez_name_agent_b}")
                        name_agent_b = entrez_name_agent_b

                # Ground agents
                agent_a = self._make_agent(name_agent_a, wormbase_id_agent_a,
                                           entrez_id_agent_a, up_id_agent_a, intact_id_agent_a) or {}
                agent_b = self._make_agent(name_agent_b, wormbase_id_agent_b,
                                           entrez_id_agent_b, up_id_agent_b, intact_id_agent_b) or {}

                # Skip any agents with no grounding
                if agent_a is None or agent_b is None:
                    continue

                # Get evidence
                pmid = None
                mint = None
                imex = None
                doi = None
                pub_id_info = self._id_conversion(wb_row.publication_identifiers) or {}
                if not pub_id_info:
                    logger.warning(f"No publication info found: {wb_row}")
                else:
                    if pub_id_info.get('pubmed'):
                        pmid = pub_id_info.get('pubmed')[0]
                    if pub_id_info.get('mint'):
                        mint = pub_id_info.get('mint')[0]
                    if pub_id_info.get('imex'):
                        imex = pub_id_info.get('imex')[0]
                    if pub_id_info.get('doi'):
                        doi = pub_id_info.get('doi')[0]

                text_refs = {}
                if pmid:
                    text_refs['PMID'] = pmid
                if doi:
                    text_refs['DOI'] = doi
                if mint:
                    text_refs['MINT'] = mint
                if imex:
                    text_refs['IMEX'] = imex

                source = None
                int_id_info = self._id_conversion(wb_row.interaction_identifiers) or {}
                if not int_id_info:
                    logger.warning(f"No interaction ID found: {wb_row}")
                else:
                    if 'wormbase' in int_id_info:
                        source = 'wormbase'
                    else:
                        key = next(iter(int_id_info), None)
                        source = (int_id_info.get(key) or [None])[0]

                source_id = (int_id_info.get(source) or [None])[0]

                # Incorporate info from the wormbase-to-entrez ID mapping file into Evidence as annotations
                full_annotations = {}
                full_annotations['interaction_info'] = wb_row._asdict()
                full_annotations['entrez_info_agent_a'] = {}
                full_annotations['entrez_info_agent_b'] = {}
                if entrez_id_agent_a:
                    full_annotations['entrez_info_agent_a'] = symbol_to_annotation_dict.get(entrez_id_agent_a) or {}
                if entrez_id_agent_b:
                    full_annotations['entrez_info_agent_b'] = symbol_to_annotation_dict.get(entrez_id_agent_b) or {}


                ev = Evidence(source_api=source,
                              source_id=source_id,
                              pmid=text_refs.get('PMID'),
                              text_refs=text_refs,
                              annotations=full_annotations
                              )
                # Make statement
                interaction_type = None
                int_type_info = self._interaction_type_conversion(wb_row.interaction_types) or {}
                if not int_type_info:
                    logger.warning(f"No interaction type found: {wb_row}")
                else:
                    if int_type_info.get('psi-mi'):
                        interaction_type = int_type_info.get('psi-mi')[0]
                    else:
                        key = next(iter(int_type_info), None)
                        interaction_type = (int_type_info.get(key) or [None])[0]

                    if 'enhancement' in interaction_type:
                        s = Activation(agent_a, agent_b, evidence=ev)
                    elif 'suppression' in interaction_type:
                        s = Inhibition(agent_a, agent_b, evidence=ev)
                    elif 'phosphorylation reaction' in interaction_type:
                        s = Phosphorylation(agent_a, agent_b, evidence=ev)
                    elif 'demethylation reaction' in interaction_type:
                        s = Demethylation(agent_a, agent_b, evidence=ev)
                    elif 'methylation reaction' in interaction_type:
                        s = Methylation(agent_a, agent_b, evidence=ev)
                    else:
                        s = Association([agent_a, agent_b], evidence=ev)

                    self.statements.append(s)

            except Exception as e:
                logger.error(f"Error occurred at row {idx} for the {interaction_type} between agents {agent_a} and "
                             f"{agent_b}: {e}")

    def _make_agent(self, symbol, wormbase_id, entrez_id, up_id, intact_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        symbol : str
            A plain text symbol, or None if not listed.
        wormbase_id : str
            WormBase identifier
        entrez_id : str
            Entrez identifier
        up_id : str
            UniProt identifier
        intact_id : str
            IntAct identifier

        Returns
        -------
        agent : indra.statements.Agent
            A grounded agent object.
        """
        db_refs = {}
        name = symbol
        if wormbase_id:
            db_refs['WB'] = wormbase_id
        if entrez_id:
            db_refs['EGID'] = entrez_id
        if up_id:
            db_refs['UP'] = up_id
        # if intact_id:
        #     db_refs['INTACT'] = intact_id
        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name

        # At the time of writing this, the name was never None but
        # just in case
        if name is None:
            return None

        return Agent(name, db_refs=db_refs)

    def _alias_conversion(self, raw_value: str):
        """Return dictionary with keys corresponding to name types and values
        to agent names by decomposing the string value in one of 'Alias(es) interactor A' or
        'Alias(es) interactor B'.

        Example string value: 'wormbase:dpy-21(public_name)|wormbase:Y59A8B.1(sequence_name)'

        Parameters
        ----------
        raw_value : str
            The raw value in 'Alias(es) interactor A' or 'Alias(es) interactor B'
            for a particular row.

        Returns
        -------
        name_info : dict
            Dictionary with name types as keys and agent names as values (for C. elegans interaction data, the
            primary name and the one used corresponds with the key 'public_name').
        """
        # import re
        if not raw_value:
            return {}

        # Remove the strings "public name" and all double quotes (only a few special cases in the data have this)
        cleaned_value = raw_value.replace('"public name: ', '').replace('"', '')
        name_info = {}
        for sub in cleaned_value.split('|'): # 'Alias(es) interactor _' can contain multiple aliases separated by "|".
            if ':' in sub and '(' in sub:
                match = re.search(r'\(([^)]+)\)', sub)  # Extract text inside parentheses
                if match:
                    key = match.group(1)
                    val = sub.split(':')[1].split('(')[0]
                    if key not in name_info:
                        name_info[key] = [val]
                    else:
                        name_info[key].append(val)
        return name_info

    def _id_conversion(self, raw_value: str):
        """Decompose the string value in columns 'ID(s) interactor A', 'ID(s) interactor B',
        'Alt. ID(s) interactor A', 'Alt. ID(s) interactor B', 'Publication ID(s)', or
        'Interaction identifier(s)' and return dictionary with keys corresponding to
        database/source names and values to identifiers.

        Example string values: 'wormbase:WBGene00006352', 'entrez gene/locuslink:178272',
        'pubmed:36969515', 'wormbase:WBInteraction000000001'.

        Parameters
        ----------
        raw_value : str
            The raw value in whichever ID column is being converted.

        Returns
        -------
        source_id_info : dict
            Dictionary with database/source names as keys and identifiers as values. Unique keys for
            'ID(s) interactor _' in C. elegans interaction data are 'wormbase' and 'entrez gene/locuslink'.
            Unique keys for 'Publication ID(s)' in C. elegans interaction data are 'pubmed'.
        """
        if not raw_value or not isinstance(raw_value, str):
            return {}
        id_info = {}
        for sub in raw_value.split('|'):
            if ':' in sub:
                key = sub.split(':')[-2]
                val = sub.split(':')[-1]
            if key not in id_info:
                id_info[key] = [val]
            else:
                id_info[key].append(val)
        return id_info

    def _interaction_type_conversion(self, raw_value: str):
        """Decompose the string value in columns 'Interaction type(s)' and return dictionary with keys
        corresponding to database/source names and values to identifiers.

        Example string values: 'wormbase:WBGene00006352', 'entrez gene/locuslink:178272',
        'pubmed:36969515', 'wormbase:WBInteraction000000001'.

        Parameters
        ----------
        raw_value : str
            The raw value in whichever ID column is being converted.

        Returns
        -------
        source_id_info : dict
            Dictionary with database/source names as keys and identifiers as values. Unique keys for
            'ID(s) interactor _' in C. elegans interaction data are 'wormbase' and 'entrez gene/locuslink'.
            Unique keys for 'Publication ID(s)' in C. elegans interaction data are 'pubmed'.
        """
        import re
        if not raw_value:
            return {}
        type_info = {}
        for sub in raw_value.split('|'):
            if all(char in sub for char in (':', '(', ')')):
                key = sub.split(':')[0]
                val = re.search(r'\((.*)\)', sub).group(1)  # Extract text inside outermost parentheses
                if key not in type_info:
                    type_info[key] = [val]
                else:
                    type_info[key].append(val)
        return type_info

    @staticmethod
    def _download_wormbase_data(url, as_df):
        """Downloads gzipped, tab-separated WormBase data in .tab2 format.

        Parameters
        -----------
        url : str
            URL of the WormBase gzip file.
        as_df : Boolean
            If True, returns a Pandas DataFrame object. If False,
            returns an iterable list of file rows. The WormBase-to-Entrez ID
            mappings file is needed in pd.DataFrame format, while the interactions
            data are needed in as an iterable list of file rows.

        Returns
        -------
        csv_reader : list
            An iterable list of rows in the file (header has already
            been skipped).
        """
        res = requests.get(url)
        if res.status_code != 200:
            raise Exception('Unable to download WormBase data: status code %s'
                            % res.status_code)

        gzip_bytes = BytesIO(res.content)
        with gzip.open(gzip_bytes, 'rt') as gz_file:
            # Locate the header line (last line that starts with '#')
            header_index = None
            for i, line in enumerate(gz_file):
                if line.startswith('#') and not line.strip().startswith('######'):
                    header_line = line.strip('#').strip()
                    header_index = i

            if header_index is None:
                raise Exception('Header not found in the file.')

            gzip_bytes.seek(0)
            gz_file = gzip.open(gzip_bytes, 'rt')

            if as_df:
                # Create a dataframe
                columns = header_line.split('\t')
                df = pd.read_csv(gz_file, sep='\t', skiprows=header_index + 1,
                                        names=columns, na_values="-")
                return df
            else:
                # Skip all rows up to and including the header
                for _ in range(header_index + 1):
                    next(gz_file)

                csv_reader = list(csv.reader(gz_file, delimiter='\t'))  # Create list of rows
                return csv_reader

    def _read_wormbase_data(self, wormbase_file=None):
        """Return a csv.reader for a TSV file.

            Returns
            -------
            csv_reader : list
                An iterable list of rows in the file (header has already
                been skipped).
            """

        file_path = wormbase_file

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            # Locate the header line (last line that starts with '#')
            with open(file_path, 'r') as file:
                header_index = None
                for i, line in enumerate(file):
                    if line.startswith('#') and not line.strip().startswith('######'):
                        header_index = i

                if header_index is None:
                    raise Exception('Header not found in the file.')

                # Skip all rows up to and including the header
                file.seek(0)
                for _ in range(header_index + 1):
                    next(file)

                csv_reader = list(csv.reader(file, delimiter='\t'))  # Create list of row
                return csv_reader

        except Exception as e:
            raise Exception(f"Error occurred while reading WormBase CSV: {e}")

