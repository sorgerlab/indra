# -*- coding: utf-8 -*-

"""An INDRA processor for miRTarBase."""

import pandas as pd

from indra.databases import hgnc_client, mirbase_client
from indra.statements import Agent, DecreaseAmount, Evidence

__all__ = [
    'MirtarbaseProcessor',
]

VERSION = '7.0'
DATA_URL = f'http://mirtarbase.mbc.nctu.edu.tw/' \
    f'cache/download/{VERSION}/miRTarBase_MTI.xlsx'


class MirtarbaseProcessor:
    """Extracts INDRA statements from miRTarBase curated interactions.

    Parameters
    ----------
    mirtarbase_file : Optional[str]
        The file path or URL to the miRTarBase file as a Excel document.
        If not provided, the miRTarBase data is downloaded from the
        miRTarBase website.

    Attributes
    ----------
    statements: list[indra.statements.DecreaseAmount]
        Extracted INDRA DecreaseAmount statements.
    """

    def __init__(self, mirtarbase_file=None):
        df = pd.read_excel(mirtarbase_file or DATA_URL)

        try:
            from tqdm import tqdm
        except ImportError:
            it = df.values
        else:
            it = tqdm(df.values, total=len(df.index))

        self.statements = []
        for row in it:
            try:
                (mirtarbase_id, mirna_name, mirna_species, gene_name,
                 entrez_id, target_species, exp, sup_type, pmid) = row
            except ValueError:
                it.write(f'Issue with row: {row}')
                continue

            try:
                entrez_id = str(int(entrez_id))
            except ValueError:
                it.write(f'Issue on {mirtarbase_id} with'
                         f' Entrez ID {entrez_id}')
                continue

            target_db_refs = self._make_target_db_refs(entrez_id)
            target_agent = Agent(gene_name, db_refs=target_db_refs)

            mirna_db_refs = self._make_mirna_db_refs(mirna_name)
            mirna_agent = Agent(mirna_name, db_refs=mirna_db_refs)

            # TODO use experimental context information
            # TODO use species information (corresponding to gene)

            evidence = Evidence(
                source_api='mirtarbase',
                source_id=mirtarbase_id,
                pmid=pmid,
            )

            statement = DecreaseAmount(
                mirna_agent,
                target_agent,
                evidence=evidence,
            )
            self.statements.append(statement)

    @staticmethod
    def _make_target_db_refs(entrez_id):
        db_refs = {
            'TEXT': entrez_id,
        }

        hgnc_id = hgnc_client.get_hgnc_from_entrez(entrez_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id
            up_id = hgnc_client.get_uniprot_id(hgnc_id)
            if up_id is not None:
                db_refs['UP'] = up_id

        return db_refs

    @staticmethod
    def _make_mirna_db_refs(mirna_name):
        db_refs = {
            'TEXT': mirna_name,
        }

        mirbase_id = \
            mirbase_client.get_mirbase_id_from_mirbase_name(mirna_name)
        if mirbase_id is not None:
            db_refs['MIRBASE'] = mirbase_id
            hgnc_id = mirbase_client.get_hgnc_id_from_mirbase_id(mirbase_id)
            if hgnc_id is not None:
                db_refs['HGNC'] = hgnc_id

        return db_refs


def main():
    mp = MirtarbaseProcessor()
    stmts = mp.statements
    print(stmts[0])

    import pickle
    import os
    path = os.path.join(os.path.expanduser('~'), 'Desktop', 'mirbase_indra.pkl')
    with open(path, 'wb') as f:
        pickle.dump(stmts, f)


if __name__ == '__main__':
    main()
