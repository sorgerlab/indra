from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['TasProcessor']

from indra.statements import Inhibition, Agent
from indra.databases.lincs_client import LincsClient


class TasProcessor(object):
    """A processor for Target Affinity Spectrum data compiled by N. Moret.

    This data was compiled in the HMS LSP as an improvement on the "arbitrary"
    selection of targets present in the similar LINCS dataset.
    """
    def __init__(self, data):
        self._data = data
        self._lc = LincsClient()

        self.statements = []
        for row in data:
            self._process_row(row)
        return

    def _process_row(self, row):
        drug = self._extract_drug(row['hms_id'])
        prot = self._extract_protein(row['gene_id'])
        ev_list = self._make_evidence(row)
        return Inhibition(drug, prot, evidence=ev_list)

    def _extract_drug(self, hms_id):
        refs = self._lc.get_small_molecule_ref(hms_id, abbreviated=True)
        name = self._lc.get_small_molecule_name(hms_id, abbreviated=True)
        return Agent(name, db_refs=refs)

    def _extract_protein(self, gene_id):
        pass

    def _make_evidence(self, row):
        pass
