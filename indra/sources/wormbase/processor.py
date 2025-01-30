## TO-DO:
## – modify the _make_agent function to handle wormbase data

import re
import csv
import os
import tqdm
import logging
import requests
from io import BytesIO, StringIO
# from zipfile import ZipFile
import gzip
from collections import namedtuple
# from indra.util import read_unicode_csv
# from indra.util import read_mitab_csv
# from indra.statements import Agent, Complex, Evidence
# from indra.ontology.standardize import standardize_name_db_refs

logger = logging.getLogger(__name__)


wormbase_file_url = ('https://fms.alliancegenome.org/download/'
                     'INTERACTION-GEN_WB.tsv.gz')


# The explanation for each column of the tsv file is here:
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

class WormBaseProcessor(object):
    """Extracts INDRA statements from WormBase interaction data.

        Parameters
        ----------
        wormbase_file : str
            The file containing the WormBase data in .tsv format. If not provided,
            the WormBase data is downloaded from the WormBase website (technically
            alliancegenome.org).
        physical_only : boolean
            If True, only physical interactions are included (e.g., genetic
            interactions are excluded). If False, all interactions are included).

        Attributes
        ----------
        statements : list[indra.statements.Statements]
            Extracted INDRA statements.
        physical_only : boolean
            Indicates whether only physical interactions were included during
            statement processing.
        """

    def __init__(self, wormbase_file=None, physical_only=True):
        self.statements = []
        self.physical_only = physical_only

        # If a path to the file is included, process it, skipping the header
        if wormbase_file:
            rows = self._read_wormbase_csv(wormbase_file)
        # If no file is provided, download from web
        else:
            logger.info('No data file specified, downloading from WormBase '
                        'at %s' % wormbase_file_url)
            rows = self._download_wormbase_data(wormbase_file_url)

        # Process the rows into Statements
        for row in tqdm.tqdm(rows, desc='Processing WormBase rows'):
            # There are some extra columns that we don't need to take and
            # thereby save space in annotations
            filt_row = [None if item == '-' else item
                        for item in row][:len(columns)]
            wb_row = _WormBaseRow(*filt_row)

            # Get names of agents using wb_row.aliases_interactor_a and aliases_interactor_b
            name_info_agent_a = self._alias_conversion(wb_row.aliases_interactor_a)
            name_agent_a = name_info_agent_a.get('public_name')
            name_info_agent_b = self._alias_conversion(wb_row.aliases_interactor_b)
            name_agent_b = name_info_agent_b.get('public_name')

            # Get db_refs using the wormbase ID and entrez ID in wb_row.ids_interactor_a and wb_row.ids_interactor_b
            db_id_info_agent_a = self._id_conversion(wb_row.ids_interactor_a)
            wormbase_id_agent_a = db_id_info_agent_a.get('wormbase')
            entrez_id_agent_a = db_id_info_agent_a.get('entrez gene/locuslink')
            db_id_info_agent_b = self._id_conversion(wb_row.ids_interactor_b)
            wormbase_id_agent_b = db_id_info_agent_b.get('wormbase')
            entrez_id_agent_b = db_id_info_agent_b.get('entrez gene/locuslink')



        #
        #     # Filter out non-physical interactions if desired
        #     # if self.physical_only and wb_row.exp_system_type != 'physical':
        #     #     continue
        #
        #     # Ground agents
        #     agent_a = self._make_agent(wb_row.symbol_a, wb_row.entrez_a,
        #                                wb_row.swissprot_a, wb_row.trembl_a)
        #     agent_b = self._make_agent(wb_row.symbol_b, wb_row.entrez_b,
        #                                wb_row.swissprot_b, wb_row.trembl_b)


    def _make_agent(self, symbol, wormbase_id, entrez_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        wormbase_id : str
            WormBase identifier
        entrez_id : str
            Entrez id number
        symbol : str
            A plain text symbol, or None if not listed.

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
        to agent names (or aliases) by decomposing the string value in Alias(es) interact A
        or Alias(es) interactor B.

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
        name_info = {}
        for sub in raw_value.split('|'): # 'Alias(es) interactor _' can contain multiple aliases separated by "|".
            if ':' in sub and '(' in sub:
                match = re.search(r'\(([^)]+)\)', sub)  # Extract text inside parentheses
                if match:
                    key = match.group(1)
                    val = sub.split(':')[1].split('(')[0]
                    name_info[key] = val
        return name_info


    def _id_conversion(self, raw_value: str):
        """Decompose the string value in columns 'ID(s) interactor A', 'ID(s) interactor B',
        or 'Publication ID(s)' and return dictionary with keys corresponding to database/source names and values
        to identifiers.

        Example string values: 'wormbase:WBGene00006352', 'entrez gene/locuslink:178272', 'pubmed:36969515'.

        Parameters
        ----------
        raw_value : str
            The raw value in 'ID(s) interactor A' or 'ID(s) interactor B'

        Returns
        -------
        source_id_info : dict
            Dictionary with database/source names as keys and identifiers as values. Unique keys for
            'ID(s) interactor _' in C. elegans interaction data are 'wormbase' and 'entrez gene/locuslink'.
            Unique keys for 'Publication ID(s)' in C. elegans interaction data are 'pubmed'.
        """
        # import re
        if not raw_value:
            return {}
        id_info = {}
        for sub in raw_value.split('|'):
            if ':' in sub:
                key = sub.split(':')[0]
                val = sub.split(':')[1]
                id_info[key] = val
        return id_info


    def _download_wormbase_data(url):
        """Downloads gzipped, tab-separated WormBase data in .tab2 format.

        Parameters
        -----------
        url : str
            URL of the WormBase gzip file.

        Returns
        -------
        csv.reader
            A csv.reader object for iterating over the rows (header has already
            been skipped).
        """
        res = requests.get(wormbase_file_url)
        if res.status_code != 200:
            raise Exception('Unable to download WormBase data: status code %s'
                            % res.status_code)

        gzip_bytes = BytesIO(res.content)
        with gzip.open(gzip_bytes, 'rt') as gz_file:  # Open the .gz file in text mode
            # Locate the header line (last line starting with '#')
            # with open(file_path, 'r') as file:
            #     lines = file.readlines()
            #
            # header_line = None
            header_index = None
            for i, line in enumerate(gz_file):
                if line.startswith('#') and not line.strip().startswith('######'):
                    # header_line = line.strip('#').strip()
                    header_index = i

            if header_index is None:
                raise Exception('Header not found in the file.')

            # Reset the file pointer to the beginning after locating the header
            gzip_bytes.seek(0)
            gz_file = gzip.open(gzip_bytes, 'rt')  # Reinitialize gz_file

            # Skip rows until the specified header_index
            for _ in range(header_index + 1):  # Skip all rows up to and including the header
                next(gz_file)

            csv_reader = csv.reader(gz_file, delimiter='\t')  # Read TSV content
            return csv_reader  # Return csv.reader for iteration

    def _read_wormbase_csv(file_path: str) -> csv.reader:
        """Return a csv.reader for a TSV file.

            Parameters
            ----------
            file_path : str
                Path to TSV file that is to be read.

            Returns
            -------
            csv_reader : csv.reader
                CSV reader for iteration.
            """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        try:
            # Locate the header line (last line starting with '#')
            with open(file_path, 'r') as file:
                header_index = None
                for i, line in enumerate(file):
                    if line.startswith('#') and not line.strip().startswith('######'):
                        # header_line = line.strip('#').strip()
                        header_index = i

                if header_index is None:
                    raise Exception('Header not found in the file.')

                # Reset the file pointer to the beginning after locating the header
                file.seek(0)

                # Skip rows until the specified header_index
                for _ in range(header_index + 1):  # Skip all rows up to and including the header
                    next(file)

                csv_reader = csv.reader(file, delimiter='\t')  # Read TSV content
                return csv_reader  # Return csv.reader for iteration
