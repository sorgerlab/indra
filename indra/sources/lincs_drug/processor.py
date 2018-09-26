from __future__ import absolute_import, print_function, unicode_literals

__all__ = ['LincsProcessor']

import re

from indra.statements import Agent, Inhibition, Evidence
from indra.databases.lincs_client import LincsClient


class LincsProcessor(object):
    """Processor for the HMS LINCS drug target dataset.

    Parameters
    ----------
    lincs_data : list[dict]
        A list of dicts with keys set by the header of the csv, and values from
        the data in the csv.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of indra statements extracted from the CSV file.
    """

    def __init__(self, lincs_data):
        self._data = lincs_data
        self._lc = LincsClient()

        # Process all the lines (skipping the header)
        self.statements = []
        for line in self._data:
            self._process_line(line)
        return

    def _process_line(self, line):
        drug = self._extract_drug(line)
        prot = self._extract_protein(line)
        evidence = self._make_evidence(line)
        self.statements.append(Inhibition(drug, prot, evidence=evidence))
        return

    def _extract_drug(self, line):
        drug_name = line['Small Molecule Name']
        lincs_id = line['Small Molecule HMS LINCS ID']
        db_refs = self._lc.get_small_molecule_ref(lincs_id)
        return Agent(drug_name, db_refs=db_refs)

    def _extract_protein(self, line):
        prot_name = line['Protein Name']
        prot_id = line['Protein HMS LINCS ID']
        sm_id = line['Small Molecule HMS LINCS ID']
        db_refs = {}
        if prot_id:
            db_refs.update(self._lc.get_protein_ref(prot_id))
        if sm_id:
            db_refs.update(self._lc.get_small_molecule_ref(sm_id))
        assert db_refs, "We didn't get any refs for: %s" % prot_name
        return Agent(prot_name, db_refs=db_refs)

    def _make_evidence(self, line):
        ev_list = []
        key_refs = line['Key References'].split(';')
        generic_notes = {
            'is_nominal': line['Is Nominal'],
            'effective_concentration': line['Effective Concentration']
            }
        patt = re.compile('(?:pmid|pubmed\s+id):\s+(\d+)', re.IGNORECASE)
        for ref in key_refs:
            # Only extracting pmids, but there is generally more info available.
            m = patt.search(ref)
            if m is None:
                pmid = None
            else:
                pmid = m.groups()[0]
            annotations = {'reference': ref}
            annotations.update(generic_notes)
            ev = Evidence('lincs', pmid=pmid, annotations=annotations,
                          epistemics={'direct': True})
            ev_list.append(ev)
        return ev_list


