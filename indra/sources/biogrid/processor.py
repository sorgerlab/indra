from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import logging
import operator
import itertools
import collections
from indra.util import read_unicode_csv
from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client

logger = logging.getLogger('biogrid')


class BiogridProcessor(object):
    """The BioGrid processor reads in a tab-seperated file with biogrid
    interactions and produces a list of INDRA statements.

    Parameters
    ----------
    filename: str
        The filename of the biogrid data file

    Attributes:
    -----------
    statements: list[indra.statements.Statement]
        Extracted INDRA statements
    """
    def __init__(self, filename):
        self.statements = []

        rows = read_unicode_csv(filename, '\t')
        for row in rows:
            # The explanation for each column of the tsv file is here:
            # https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_2.0
            source_id = row[0]
            entrez_a = row[1]
            entrez_b = row[2]
            text_a = row[5]
            text_b = row[6]
            pmid = row[14]
            system_type = row[12]

            if system_type != 'physical':
                continue

            # Ground agents
            agent_a = self._make_agent(entrez_a, text_a)
            agent_b = self._make_agent(entrez_b, text_b)

            # Evidence
            ev = Evidence(source_api='biogrid',
                          source_id=source_id,
                          pmid=pmid,
                          text=None,
                          annotations={'tsv_row': row})

            # Make statement
            s = Complex([agent_a, agent_b], evidence=ev)
            self.statements.append(s)

    def _make_agent(self, entrez_id, text_id):
        """Make an Agent object, appropriately grounded.

        Parameters
        ----------
        entrez_id: str
            Entrez id number
        text_id: str
            A plain text systematic name, or - if not listed

        Returns
        -------
        agent: indra.statements.Agent
            A grounded agent object
        """
        hgnc_name, db_refs = self._make_db_refs(entrez_id, text_id)
        if hgnc_name is not None:
            name = hgnc_name
        else:
            name = text_id

        return Agent(name, db_refs=db_refs)

    def _make_db_refs(self, entrez_id, text_id):
        """Looks up the hgnc id and name, as well as the uniprot id.

        Parameters
        ----------
        entrez_id: str
            Entrez id number
        text_id:
            A plain text systematic name, or - if not listed in biogrid

        Returns
        -------
        hgnc_name: str
            name from the HGNC database
        db_refs: dict
            db_refs grounding dictionary, used when constructing the Agent
            object
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
