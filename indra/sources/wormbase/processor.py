## TO-DO:
## – modify the _make_agent function so that standardized

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
from indra.statements import Agent, Evidence, Activation, Association, Inhibition
from indra.ontology.standardize import standardize_name_db_refs

# comment

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

    def __init__(self, wormbase_file=None):
        self.statements = []
        self.wormbase_file = wormbase_file

        # If a path to the file is included, process it, skipping the header
        if self.wormbase_file:
            rows = self._read_wormbase_data()
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

            # Get db_refs using the wormbase ID and entrez ID in wb_row.ids_interactor_a
            # and wb_row.ids_interactor_b
            db_id_info_agent_a = self._id_conversion(wb_row.ids_interactor_a)
            wormbase_id_agent_a = db_id_info_agent_a.get('wormbase')
            entrez_id_agent_a = db_id_info_agent_a.get('entrez gene/locuslink')
            db_id_info_agent_b = self._id_conversion(wb_row.ids_interactor_b)
            wormbase_id_agent_b = db_id_info_agent_b.get('wormbase')
            entrez_id_agent_b = db_id_info_agent_b.get('entrez gene/locuslink')

            # Ground agents
            agent_a = self._make_agent(name_agent_a, wormbase_id_agent_a, entrez_id_agent_a)
            agent_b = self._make_agent(name_agent_b, wormbase_id_agent_b, entrez_id_agent_b)

            # Skip any agents with neither HGNC grounding or string name
            if agent_a is None or agent_b is None:
                continue
            # Get evidence
            pmid = self._id_conversion(wb_row.publication_identifiers).get('pubmed')
            doi = self._id_conversion(wb_row.publication_identifiers).get('doi')
            text_refs = {}
            if pmid:
                text_refs['PMID'] = pmid
            elif doi:
                text_refs['DOI'] = doi

            source_id = self._id_conversion(wb_row.interaction_identifiers).get('wormbase')
            interaction_annotations = self._id_conversion(wb_row.interaction_annotations).get('wormbase')
            ev = Evidence(source_api='wormbase',
                          source_id=source_id,
                          pmid=text_refs.get('PMID'),
                          text_refs=text_refs,
                          # annotations=interaction_annotations
                          annotations=dict(wb_row._asdict())
                          )
            # Make statement
            interaction_type = self._interaction_type_conversion(wb_row.interaction_types).get('psi-mi')
            if 'enhancement' in interaction_type:
                s = Activation([agent_a, agent_b], evidence=ev)
            elif 'suppression' in interaction_type:
                s = Inhibition([agent_a, agent_b], evidence=ev)
            else:
                s = Association([agent_a, agent_b], evidence=ev)
            self.statements.append(s)

    def _make_agent(self, symbol, wormbase_id, entrez_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        symbol : str
            A plain text symbol, or None if not listed.
        wormbase_id : str
            WormBase identifier
        entrez_id : str
            Entrez id number

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
        #
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
        'Publication ID(s)', or 'Interaction identifier(s)' and return dictionary with keys
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
        if not raw_value:
            return {}
        id_info = {}
        for sub in raw_value.split('|'):
            if ':' in sub:
                key = sub.split(':')[0]
                val = sub.split(':')[1]
                id_info[key] = val
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
                type_info[key] = val
        return type_info

    @staticmethod
    def _download_wormbase_data(url):
        """Downloads gzipped, tab-separated WormBase data in .tab2 format.

        Parameters
        -----------
        url : str
            URL of the WormBase gzip file.

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
                    header_index = i

            if header_index is None:
                raise Exception('Header not found in the file.')

            gzip_bytes.seek(0)
            gz_file = gzip.open(gzip_bytes, 'rt')

            # Skip all rows up to and including the header
            for _ in range(header_index + 1):
                next(gz_file)

            csv_reader = list(csv.reader(gz_file, delimiter='\t'))  # Create list of rows
            return csv_reader

    def _read_wormbase_data(self):
        """Return a csv.reader for a TSV file.

            Returns
            -------
            csv_reader : list
                An iterable list of rows in the file (header has already
                been skipped).
            """

        file_path = self.wormbase_file

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

