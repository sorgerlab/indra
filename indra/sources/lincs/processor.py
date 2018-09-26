from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['LincsProcessor']

import re

from indra.statements import Agent, Inhibition, Evidence
from indra.sources.lincs.lincs_client import get_small_molecule_data, \
    get_protein_data


class LincsProcessor(object):
    """Processor for the HMS LINCS drug target dataset.

    Parameters
    ----------
    lincs_csv : str
        A csv of the lincs data.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of indra statements extracted from the CSV file.
    """

    def __init__(self, lincs_data):
        self._data = lincs_data
        self._sm_data = get_small_molecule_data()
        self._prot_data = get_protein_data()

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
        db_refs = _build_db_refs(lincs_id, self._sm_data[lincs_id],
                                 chembl='ChEMBL ID', chebi='ChEBI ID',
                                 pubchem='PubChem CID', inchi='InChi Key',
                                 lincs='LINCS ID')
        return Agent(drug_name, db_refs=db_refs)

    def _extract_protein(self, line):
        prot_name = line['Protein Name']
        prot_id = line['Protein HMS LINCS ID']
        sm_id = line['Small Molecule HMS LINCS ID']
        db_refs = {}
        if prot_id:
            # TODO: We could get phosphorylation states from the prtein data.
            db_refs.update(_build_db_refs(prot_id, self._prot_data[prot_id],
                                          hgnc='Gene ID', uniprot='UniProt ID'))
        if sm_id:
            db_refs.update(_build_db_refs(sm_id, self._sm_data[sm_id],
                                          chembl='ChEMBL ID', chebi='ChEBI ID',
                                          pubchem='PubChem CID',
                                          inchi='InChi Key', lincs='LINCS ID'))
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
            ev = Evidence('lincs', pmid=pmid, annotations=annotations)
            ev_list.append(ev)
        return ev_list


def _build_db_refs(lincs_id, data, **mappings):
    db_refs = {'HMS-LINCS': lincs_id}
    for db_ref, key in mappings.items():
        if data[key]:
            db_refs[db_ref.upper()] = data[key]
    return db_refs
