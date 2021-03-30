# -*- coding: utf-8 -*-

"""Process statements from the `Drug Gene Interaction Database (DGI-DB) <http://www.dgidb.org>`_."""

import logging
from typing import Iterable, List, Optional, Type

import pandas as pd

import pystow
from indra.statements import Activation, Agent, DecreaseAmount, Evidence, IncreaseAmount, Inhibition, Statement

__all__ = [
    'DGIProcessor',
]

logger = logging.getLogger(__name__)


class DGIProcessor:
    """Processor to extract INDRA Statements from DGI content.

    Parameters
    ----------
    version : str
        The version of DGI to use
    """

    #: A list of INDRA Statements that were extracted from DGI content.
    statements: List[Statement]

    def __init__(self, version: Optional[str] = None):
        self.df = get_df(version=version)
        self.statements = []

    def extract_statements(self) -> List[Statement]:
        for gene_name, ncbigene_id, interaction, drug_name, drug_curie, pmids in self.df.values:
            statements = self.row_to_statements(gene_name, ncbigene_id, interaction, drug_name, drug_curie, pmids)
            self.statements.extend(statements)
        return self.statements

    @classmethod
    def row_to_statements(
        cls,
        gene_name,
        ncbigene_id,
        interactions,
        drug_name,
        drug_curie,
        pmids,
    ) -> Iterable[Statement]:
        gene_agent = Agent.from_refs(gene_name, {'EGID': ncbigene_id})

        try:
            drug_namespace, drug_identifier = drug_curie.split(':', 1)
            drug_namespace = drug_namespace.upper()
        except ValueError:
            logger.warning('could not parse drug CURIE: %s', drug_curie)
            return

        drug_agent = Agent.from_refs(drug_name, {drug_namespace: drug_identifier})

        evidence = [
            Evidence(source_api='dgi', pmid=pmid, annotations={
                'interaction': interactions,
            })
            for pmid in pmids or [None]
        ]
        statement_cls = _get_statement_type(interactions)
        if statement_cls is not None:
            yield statement_cls(drug_agent, gene_agent, evidence=evidence)


def get_df(version: Optional[str] = None) -> pd.DataFrame:
    """Get the DGI interaction dataframe."""
    if version is None:
        version = '2021-Jan'  # TODO use bioversions to look up
    url = f'https://www.dgidb.org/data/monthly_tsvs/{version}/interactions.tsv'
    df = pystow.ensure_csv(
        'indra', 'sources', 'dgi', version,
        url=url,
        read_csv_kwargs=dict(
            usecols=['gene_name', 'entrez_id', 'interaction_types', 'drug_name', 'drug_concept_id', 'PMIDs'],
            dtype=str,
        ),
    )
    # remove rows with missing information
    df = df[df['entrez_id'].notna()]
    df = df[df['interaction_types'].notna()]
    df = df[df['drug_concept_id'].notna()]
    df['PMIDs'] = df['PMIDs'].apply(_safe_split)
    return df


def _safe_split(s: str) -> List[str]:
    if pd.isna(s):
        return []
    return [x.strip() for x in s.split(',')]


ACTIVATES_TYPES = {
    'positive allosteric modulator',
    'positive modulator',
    'activator',
    'stimulator',
    'inducer',
    'cofactor',
    'antagonist,inducer',
    'agonist,stimulator',
    'agonist,activator',
    'potentiator,activator',
    'activator,inducer',
    'blocker,activator',
    'activator,channel blocker',
    'activator,antagonist',
    'agonist,inducer',
    'agonist,potentiator',
    'modulator,activator',
    'modulator,cofactor',
    'agonist,positive modulator',
    'antagonist,potentiator',
    'modulator,inducer',
    'binder,activator',
    'inducer,substrate',
    'ligand,inducer',
    'potentiator,binder',
}

INHIBITS_TYPES = {
    'inhibitor',
    'ligand,inhibitor',
    'antibody,inhibitor',
    'inhibitor,antibody',
    'stimulator,inhibitor',
    'negative modulator,agonist,antagonist',
    'inhibitor,substrate',
    'agonist,inhibitor',
    'blocker',
    'binder,inhibitor',
    'channel blocker',
    'antagonist,inhibitor',
    'blocker,inhibitor',
    'suppressor',
    'gating inhibitor',
    'channel blocker,gating inhibitor',
    'inhibitory allosteric modulator',
    'inhibitory allosteric modulator,antagonist',
    'negative modulator',
    'negative modulator,inhibitor',
    'negative modulator,antagonist',
    'allosteric modulator,antagonist',
    'antagonist,allosteric modulator',
    'negative modulator,agonist',
    'negative modulator,inhibitor,binder',
    'antagonist,blocker',
    'modulator,inhibitor',
    'negative modulator,agonist,inhibitor',
}
INCREASE_AMOUNT_TYPES = {
    'chaperone',
}
DECREASE_AMOUNT_TYPES = {
    'antisense',
    'antisense oligonucleotide',
    'cleavage',
}
REGULATES_TYPES = {
    'agonist',  # while an agonist does the same as the native ligand, it is not inherently activate or inhibit
    'partial agonist',
    'agonist,partial agonist',
    'antagonist',  # while an agonist does the opposite as the native ligand, it is not inherently activate or inhibit
    'antagonist,partial agonist',
    'partial antagonist',
    'agonist,modulator',
    'antagonist,ligand,partial agonist',
    'agonist,allosteric modulator',
    'inverse agonist',
    'antibody',
    'modulator',
    'antagonist,binder',
    'allosteric modulator',
    'agonist,antagonist',
    'ligand',
    'antagonist,ligand',
    'modulator,antagonist',
    'antagonist,inverse agonist',
    'antagonist,multitarget',
    'antagonist,substrate',
    'modulator,ligand',
    'antagonist,antibody',
    # Contradictions
    'inhibitor,activator',
    'potentiator,inhibitor',
    'inhibitor,inducer',
}
BINDS_TYPES = {
    'binder',
    'adduct',
    'substrate',

}
SKIP_TYPES = {
    'product of',
    'product of,substrate',
    'multitarget',
    'vaccine',
}

_UNHANDLED = set()


def _get_statement_type(s: str) -> Optional[Type[Statement]]:
    if s in ACTIVATES_TYPES:
        return Activation
    if s in INCREASE_AMOUNT_TYPES:
        return IncreaseAmount
    if s in INHIBITS_TYPES:
        return Inhibition
    if s in DECREASE_AMOUNT_TYPES:
        return DecreaseAmount
    if s in REGULATES_TYPES:
        return
    if s in BINDS_TYPES or s in SKIP_TYPES:
        return
    if s not in _UNHANDLED:
        _UNHANDLED.add(s)
        logger.warning('unhandled interaction type: %s', s)


def main():
    processor = DGIProcessor()
    statements = processor.extract_statements()
    print(f'got {len(statements)} statements. Example:')
    print(statements[0])


if __name__ == '__main__':
    main()
