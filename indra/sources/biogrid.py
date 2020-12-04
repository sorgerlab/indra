import re
import csv
import tqdm
import logging
import requests
from io import BytesIO, StringIO
from zipfile import ZipFile
from collections import namedtuple
from indra.util import read_unicode_csv
from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
from indra.ontology.standardize import standardize_name_db_refs

logger = logging.getLogger(__name__)


biogrid_file_url = ('https://downloads.thebiogrid.org/Download/BioGRID/'
                    'Latest-Release/BIOGRID-ALL-LATEST.tab3.zip')


# The explanation for each column of the tsv file is here:
# https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_3.0
columns = ['biogrid_int_id',
           'entrez_a', 'entrez_b',
           'biogrid_a', 'biogrid_b',
           'syst_name_a', 'syst_name_b',
           'symbol_a', 'symbol_b',
           'syn_a', 'syn_b',
           'exp_system', 'exp_system_type',
           'author', 'publication',
           'organism_a', 'organism_b',
           'throughput', 'score', 'modification',
           'qualifications', 'tags', 'source_db',
           'swissprot_a', 'trembl_a', 'refseq_a',
           'swissprot_b', 'trembl_b', 'refseq_b']
_BiogridRow = namedtuple('BiogridRow', columns)


class BiogridProcessor(object):
    """Extracts INDRA Complex statements from Biogrid interaction data.

    Parameters
    ----------
    biogrid_file : str
        The file containing the Biogrid data in .tab2 format. If not provided,
        the BioGrid data is downloaded from the BioGrid website.
    physical_only : boolean
        If True, only physical interactions are included (e.g., genetic
        interactions are excluded). If False, all interactions are included).

    Attributes
    ----------
    statements : list[indra.statements.Statements]
        Extracted INDRA Complex statements.
    physical_only : boolean
        Indicates whether only physical interactions were included during
        statement processing.
    """
    def __init__(self, biogrid_file=None, physical_only=True):
        self.statements = []
        self.physical_only = physical_only

        # If a path to the file is included, process it, skipping the header
        if biogrid_file:
            rows = read_unicode_csv(biogrid_file, '\t', skiprows=1)
        # If no file is provided, download from web
        else:
            logger.info('No data file specified, downloading from BioGrid '
                        'at %s' % biogrid_file_url)
            rows = _download_biogrid_data(biogrid_file_url)

        # Process the rows into Statements
        for row in tqdm.tqdm(rows, desc='Processing BioGRID rows'):
            # There are some extra columns that we don't need to take and
            # thereby save space in annotations
            filt_row = [None if item == '-' else item
                        for item in row][:len(columns)]
            bg_row = _BiogridRow(*filt_row)
            # Filter out non-physical interactions if desired
            if self.physical_only and bg_row.exp_system_type != 'physical':
                continue
            # Ground agents
            agent_a = self._make_agent(bg_row.symbol_a, bg_row.entrez_a,
                                       bg_row.swissprot_a, bg_row.trembl_a)
            agent_b = self._make_agent(bg_row.symbol_b, bg_row.entrez_b,
                                       bg_row.swissprot_b, bg_row.trembl_b)
            # Skip any agents with neither HGNC grounding or string name
            if agent_a is None or agent_b is None:
                continue
            # Get evidence
            pmid_match = re.match(r'PUBMED:(\d+)',
                                  bg_row.publication)
            doi_match = re.match(r'DOI:(.*)', bg_row.publication)
            text_refs = {}
            if pmid_match:
                text_refs['PMID'] = pmid_match.groups()[0]
            elif doi_match:
                text_refs['DOI'] = doi_match.groups()[0]
            ev = Evidence(source_api='biogrid',
                          source_id=bg_row.biogrid_int_id,
                          pmid=text_refs.get('PMID'),
                          text_refs=text_refs,
                          annotations=dict(bg_row._asdict()))
            # Make statement
            s = Complex([agent_a, agent_b], evidence=ev)
            self.statements.append(s)

    def _make_agent(self, symbol, entrez_id, swissprot_id, trembl_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        entrez_id : str
            Entrez id number
        swissprot_id : str
            Swissprot (reviewed UniProt) ID.
        trembl_id : str
            Trembl (unreviewed UniProt) ID.
        symbol : str
            A plain text symbol, or None if not listed.

        Returns
        -------
        agent : indra.statements.Agent
            A grounded agent object.
        """
        db_refs = {}
        name = symbol
        if swissprot_id:
            if '|' not in swissprot_id:
                db_refs['UP'] = swissprot_id
        elif trembl_id:
            if '|' not in trembl_id:
                db_refs['UP'] = trembl_id

        if entrez_id:
            db_refs['EGID'] = entrez_id
            if 'UP' not in db_refs:
                hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
                if hgnc_id:
                    db_refs['HGNC'] = hgnc_id
        standard_name, db_refs = standardize_name_db_refs(db_refs)
        if standard_name:
            name = standard_name

        # At the time of writing this, the name was never None but
        # just in case
        if name is None:
            return None

        return Agent(name, db_refs=db_refs)


def _download_biogrid_data(url):
    """Downloads zipped, tab-separated Biogrid data in .tab2 format.

    Parameters:
    -----------
    url : str
        URL of the BioGrid zip file.

    Returns
    -------
    csv.reader
        A csv.reader object for iterating over the rows (header has already
        been skipped).
    """
    res = requests.get(biogrid_file_url)
    if res.status_code != 200:
        raise Exception('Unable to download Biogrid data: status code %s'
                        % res.status_code)
    zip_bytes = BytesIO(res.content)
    zip_file = ZipFile(zip_bytes)
    zip_info_list = zip_file.infolist()
    # There should be only one file in this zip archive
    if len(zip_info_list) != 1:
        raise Exception('There should be exactly zipfile in BioGrid zip '
                        'archive: %s' % str(zip_info_list))
    unzipped_bytes = zip_file.read(zip_info_list[0]) # Unzip the file
    biogrid_str = StringIO(unzipped_bytes.decode('utf8')) # Make file-like obj
    csv_reader = csv.reader(biogrid_str, delimiter='\t') # Get csv reader
    next(csv_reader) # Skip the header
    return csv_reader
