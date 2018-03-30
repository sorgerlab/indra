from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import csv
import logging
import itertools
import requests
from io import BytesIO, StringIO
from zipfile import ZipFile
from collections import namedtuple
from indra.util import read_unicode_csv
from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client

logger = logging.getLogger('biogrid')


biogrid_file_url = 'https://downloads.thebiogrid.org/Download/BioGRID/' + \
        'Release-Archive/BIOGRID-3.4.158/BIOGRID-ALL-3.4.158.tab2.zip'


# The explanation for each column of the tsv file is here:
# https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_2.0
_BiogridRow = namedtuple('BiogridRow',
                        ['biogrid_int_id',
                         'entrez_a', 'entrez_b',
                         'biogrid_a', 'biogrid_b',
                         'syst_name_a', 'syst_name_b',
                         'hgnc_a', 'hgnc_b',
                         'syn_a', 'syn_b',
                         'exp_system', 'exp_system_type',
                         'author', 'pmid',
                         'organism_a', 'organism_b',
                         'throughput', 'score', 'modification',
                         'phenotypes', 'qualifications', 'tags',
                         'source_db'])


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
    statements: list[indra.statements.Statements]
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
        for row in rows:
            filt_row = [None if item == '-' else item for item in row]
            bg_row = _BiogridRow(*filt_row)
            # Filter out non-physical interactions if desired
            if self.physical_only and bg_row.exp_system_type != 'physical':
                continue
            # Ground agents
            agent_a = self._make_agent(bg_row.entrez_a, bg_row.syst_name_a)
            agent_b = self._make_agent(bg_row.entrez_b, bg_row.syst_name_b)
            # Skip any agents with neither HGNC grounding or string name
            if agent_a is None or agent_b is None:
                continue
            # Get evidence
            ev = Evidence(source_api='biogrid',
                          source_id=bg_row.biogrid_int_id,
                          pmid=bg_row.pmid,
                          text=None,
                          annotations=dict(bg_row._asdict()))
            # Make statement
            s = Complex([agent_a, agent_b], evidence=ev)
            self.statements.append(s)

    def _make_agent(self, entrez_id, text_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        entrez_id : str
            Entrez id number
        text_id : str
            A plain text systematic name, or None if not listed.

        Returns
        -------
        agent : indra.statements.Agent
            A grounded agent object.
        """
        hgnc_name, db_refs = self._make_db_refs(entrez_id, text_id)
        if hgnc_name is not None:
            name = hgnc_name
        elif text_id is not None:
            name = text_id
        # Handle case where the name is None
        else:
            return None

        return Agent(name, db_refs=db_refs)

    def _make_db_refs(self, entrez_id, text_id):
        """Looks up the HGNC ID  and name, as well as the Uniprot ID.

        Parameters
        ----------
        entrez_id : str
            Entrez gene ID.
        text_id : str or None
            A plain text systematic name, or None if not listed in the
            Biogrid data.

        Returns
        -------
        hgnc_name : str
            Official HGNC symbol for the gene.
        db_refs : dict
            db_refs grounding dictionary, used when constructing the Agent
            object.
        """
        db_refs = {}
        if text_id != '-':
            db_refs['TEXT'] = text_id

        hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            if up_id is not None:
                db_refs['UP'] = up_id
        return (hgnc_name, db_refs)


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


