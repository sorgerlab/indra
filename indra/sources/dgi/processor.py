# -*- coding: utf-8 -*-

"""Processor for the `Drug Gene Interaction DB<http://www.dgidb.org>`_."""

import logging
from typing import Iterable, List, Optional, Set, Type

import pandas as pd

from ...ontology.standardize import get_standard_agent
from ...statements import (
    Activation,
    Complex,
    DecreaseAmount,
    Evidence,
    IncreaseAmount,
    Inhibition,
    Statement,
)

__all__ = [
    "DGIProcessor",
]

logger = logging.getLogger(__name__)


class DGIProcessor:
    """Processor to extract INDRA Statements from DGI content.

    Parameters
    ----------
    df : pd.DataFrame
        A pandas DataFrame for the DGI interactions file. If none given, the
        most recent version will be automatically looked up.
    version : str
        The optional version of DGI to use. If no ``df`` is given, this is also
        automatically looked up.
    """

    #: A list of INDRA Statements that were extracted from DGI content.
    statements: List[Statement]

    def __init__(
        self,
        df: Optional[pd.DataFrame] = None,
        version: Optional[str] = None,
        skip_databases: Optional[Set[str]] = None,
    ):
        if df is None:
            from .api import get_version_df
            self.version, df = get_version_df(version=version)
        else:
            self.version = version
        self.df = process_df(df)
        self.statements = []
        self.skip_databases = (
            {'DrugBank'}
            if skip_databases is None else
            set(skip_databases)
        )
        self.skipped = 0

    def extract_statements(self) -> List[Statement]:
        """Extract statements from DGI."""
        for (
            gene_name,
            ncbigene_id,
            source,
            interaction,
            drug_name,
            drug_curie,
            pmids,
        ) in self.df.values:
            if source in self.skip_databases:
                continue
            self.statements.extend(self.row_to_statements(
                gene_name,
                ncbigene_id,
                source,
                interaction,
                drug_name,
                drug_curie,
                pmids,
            ))
        return self.statements

    def row_to_statements(
        self,
        gene_name,
        ncbigene_id,
        source,
        interactions,
        drug_name,
        drug_curie,
        pmids,
    ) -> Iterable[Statement]:
        """Convert a row in the DGI dataframe to a statement."""
        gene_agent = get_standard_agent(gene_name, {"EGID": ncbigene_id})

        try:
            drug_namespace, drug_identifier = drug_curie.split(":", 1)
            drug_namespace = drug_namespace.upper()
        except ValueError:
            logger.warning("could not parse drug CURIE: %s", drug_curie)
            return

        drug_agent = get_standard_agent(
            drug_name, {drug_namespace: drug_identifier},
        )

        annotations = {
            "interactions": interactions,
            "source": source,
        }
        if self.version:
            annotations["version"] = self.version

        evidence = [
            Evidence(source_api="dgi", pmid=pmid, annotations=annotations)
            for pmid in pmids or [None]
        ]
        statement_cls = _get_statement_type(interactions)
        if statement_cls is not None:
            yield statement_cls(drug_agent, gene_agent, evidence=evidence)
        else:
            self.skipped += 1


def process_df(df: pd.DataFrame) -> pd.DataFrame:
    """Process the DGI interactions dataframe."""
    # remove rows with missing information
    df = df[df["entrez_id"].notna()]
    df = df[df["drug_concept_id"].notna()]
    df["PMIDs"] = df["PMIDs"].map(_safe_split)
    return df


def _safe_split(s: str) -> List[str]:
    if not s or pd.isna(s):
        return []
    return [x.strip() for x in s.split(",")]


ACTIVATES_TYPES = {
    "positive allosteric modulator",
    "positive modulator",
    "activator",
    "stimulator",
    "inducer",
    "cofactor",
    "antagonist,inducer",
    "agonist,stimulator",
    "agonist,activator",
    "potentiator",
    "potentiator,activator",
    "activator,inducer",
    "blocker,activator",
    "activator,channel blocker",
    "activator,antagonist",
    "agonist,inducer",
    "agonist,potentiator",
    "modulator,activator",
    "modulator,cofactor",
    "agonist,positive modulator",
    "antagonist,potentiator",
    "modulator,inducer",
    "binder,activator",
    "inducer,substrate",
    "ligand,inducer",
    "potentiator,binder",
}

INHIBITS_TYPES = {
    "inhibitor",
    "ligand,inhibitor",
    "antibody,inhibitor",
    "inhibitor,antibody",
    "stimulator,inhibitor",
    "negative modulator,agonist,antagonist",
    "inhibitor,substrate",
    "agonist,inhibitor",
    "blocker",
    "binder,inhibitor",
    "channel blocker",
    "antagonist,inhibitor",
    "blocker,inhibitor",
    "suppressor",
    "gating inhibitor",
    "channel blocker,gating inhibitor",
    "inhibitory allosteric modulator",
    "inhibitory allosteric modulator,antagonist",
    "negative modulator",
    "negative modulator,inhibitor",
    "negative modulator,antagonist",
    "allosteric modulator,antagonist",
    "antagonist,allosteric modulator",
    "negative modulator,agonist",
    "negative modulator,inhibitor,binder",
    "antagonist,blocker",
    "modulator,inhibitor",
    "negative modulator,agonist,inhibitor",
}
INCREASE_AMOUNT_TYPES = {
    "chaperone",
}
DECREASE_AMOUNT_TYPES = {
    "antisense",
    "antisense oligonucleotide",
    "cleavage",
}
REGULATES_TYPES = {
    # while an agonist does the same as the native ligand,
    # it is not inherently activate or inhibit
    "agonist",
    "partial agonist",
    "agonist,partial agonist",
    # while an agonist does the opposite as the native ligand,
    # it is not inherently activate or inhibit
    "antagonist",
    "antagonist,partial agonist",
    "partial antagonist",
    "agonist,modulator",
    "antagonist,ligand,partial agonist",
    "agonist,allosteric modulator",
    "inverse agonist",
    "antibody",
    "modulator",
    "antagonist,binder",
    "allosteric modulator",
    "agonist,antagonist",
    "antagonist,ligand",
    "modulator,antagonist",
    "antagonist,inverse agonist",
    "antagonist,multitarget",
    "antagonist,substrate",
    "modulator,ligand",
    "antagonist,antibody",
    # Contradictions
    "inhibitor,activator",
    "potentiator,inhibitor",
    "inhibitor,inducer",
    "antagonist,agonist",
}
BINDS_TYPES = {
    "ligand",
    "binder",
    "adduct",
    "substrate",
}
SKIP_TYPES = {
    "product of",
    "product of,substrate",
    "multitarget",
    "vaccine",
    "nan",
}

_UNHANDLED = set()


def _complex(a, b, evidence):
    return Complex([a, b], evidence=evidence)


def _get_statement_type(s: str) -> Optional[Type[Statement]]:
    if s in SKIP_TYPES:
        return
    if s in ACTIVATES_TYPES:
        return Activation
    if s in INCREASE_AMOUNT_TYPES:
        return IncreaseAmount
    if s in INHIBITS_TYPES:
        return Inhibition
    if s in DECREASE_AMOUNT_TYPES:
        return DecreaseAmount
    if s in REGULATES_TYPES:
        return _complex
    if s in BINDS_TYPES:
        return _complex
    if s not in _UNHANDLED:
        _UNHANDLED.add(s)
        logger.warning("unhandled interaction type: %s", s)
