"""Processors for CREEDS data."""

from copy import copy
from typing import Any, Iterable, List, Literal, Mapping, Optional, Tuple, Type

import pandas as pd
import pystow
from tqdm import tqdm

from indra import statements
from indra.databases import hgnc_client
from indra.ontology.bio import bio_ontology
from indra.ontology.standardize import get_standard_agent
from indra.sources.utils import Processor
from indra.statements import Evidence, Statement

__all__ = [
    "GeneProcessor",
]

CREEDS_MODULE = pystow.module("bio", "creeds")
BASE_URL = "http://amp.pharm.mssm.edu/CREEDS/download"
GENE_PERTURBATIONS_METADATA_URL = f'{BASE_URL}/single_gene_perturbations-v1.0.csv'
GENE_PERTURBATIONS_DATA_URL = f'{BASE_URL}/single_gene_perturbations-v1.0.json'

#: Organism label to NCBI Taxonomy Identifier
ORGANISMS = {
    "mouse": "10090",
    "human": "9606",
    "rat": "10116",
}
#: A mapping of strings used in the "pert_type" entries in
#: CREEDS data to normalized keys. Several are not curated
#: because they do not readily map to an INDRA statement
PERTURBATIONS = {
    "ko": "knockout",
    "deletion": "knockout",
    "kd": "knockdown",
    "null mutation": "knockout",
    "deficiency (mutation)": "knockdown",
    "silencing": "knockdown",
    "heterozygotic knockout (nf1+/-)": "knockout",

    "induction": "increase",
    "knock-in": "increase",
    "oe": "increase",
    "overexpression": "increase",
    "stimulation of gene product": "increase",

    "agonist activation": "activation",
    "drugactivation": "activation",
    "activemutant": "activation",
    "activation (deltanb-cateniner transgenics)": "activation",

    "druginhibition": "inhibition",
    "inhibition": "inhibition",
    "inactivation  (ikk inhibition)": "inhibition",
    "deficiency": "inhibition",
    "depletion": "inhibition",
    "depletion - sirna simut12": "inhibition",
    "small molecule inhibition": "inhibition",
    "defectivemutant": "inhibition",
    "inactivation": "inhibition",
    "drug": "inhibition",

    # misc
    "mutant": "mutation",
    "rd1 mutation": "mutation",
    "natural variation": "mutation",
    "mutation (htt exon1  142q)": "mutation",
    "g93a mutation": "mutation",
}

#: A mapping of perturbation types to statement types for the target
#: genes whose gene expression increases from the given perturbation
#: of the subject gene
UP_MAP: Mapping[str, Type[statements.RegulateAmount]] = {
    'knockout': statements.DecreaseAmount,
    'knockdown': statements.DecreaseAmount,
    'inhibition': statements.DecreaseAmount,
    'increase': statements.IncreaseAmount,
    'activation': statements.IncreaseAmount,
}
#: A mapping of perturbation types to statement types for the target
#: genes whose gene expression decreases from the given perturbation
#: of the subject gene
DOWN_MAP: Mapping[str, Type[statements.RegulateAmount]] = {
    'knockout': statements.IncreaseAmount,
    'knockdown': statements.IncreaseAmount,
    'inhibition': statements.IncreaseAmount,
    'increase': statements.DecreaseAmount,
    'activation': statements.DecreaseAmount,
}


def _process_pert_type(s: str) -> str:
    x = s.strip().lower()
    return PERTURBATIONS.get(x, x)


def preprocess_metadata(force: bool = False) -> pd.DataFrame:
    # TODO not technically necessary, but good for reference.
    #  Might delete before finishing the PR
    metadata = CREEDS_MODULE.ensure_csv(url=GENE_PERTURBATIONS_METADATA_URL, force=force, read_csv_kwargs=dict(sep=","))
    metadata = metadata[metadata.pert_type.notna()]
    metadata.id = metadata.id.map(lambda s: s[len("gene:"):])
    metadata.organism = metadata.id.map(ORGANISMS, na_action="ignore")
    metadata.pert_type = metadata.pert_type.map(_process_pert_type, na_action="ignore")
    return metadata


def _get_genes(
    record: Mapping[str, Any],
    prefix: str,
    key: Literal["down_genes", "up_genes"],
) -> List[Tuple[str, str, str]]:
    rv = []
    #: A list of 2-tuples with the gene symbol then the expression value
    expressions = record[key]
    for symbol, _ in expressions:
        try:
            _, identifier = bio_ontology.get_id_from_name(prefix, symbol)
        except TypeError:
            # "cannot unpack non-iterable NoneType object"
            # This happens if none is returned and it can't unpack the tuple
            continue
        else:
            rv.append((prefix, identifier, symbol))
    return rv


def _rat_name_from_human_name(human_name: str) -> Optional[str]:
    _, human_id = bio_ontology.get_id_from_name("HGNC", human_name)
    rat_id = hgnc_client.get_rat_id(human_id)
    return bio_ontology.get_name("RGD", rat_id)


class SimpleProcesor(Processor):
    statements: List[Statement]

    def __init__(self):
        self.statements = []

    def iter_statements(self) -> Iterable[Statement]:
        raise NotImplementedError

    def extract_statements(self) -> List[Statement]:
        if not self.statements:
            self.statements = list(self.iter_statements())
        return self.statements


class GeneProcessor(SimpleProcesor):
    """A processor for single gene perturbation experiments in CREEDS."""

    name = "creeds"

    def iter_statements(self) -> Iterable[Statement]:
        records = CREEDS_MODULE.ensure_json(url=GENE_PERTURBATIONS_DATA_URL)
        for record in tqdm(records, desc='Processing CREEDS'):
            yield from self.process_record(record)

    @staticmethod
    def process_record(record: Mapping[str, Any]) -> Iterable[Statement]:
        organism = record["organism"]
        perturbed_ncbigene_id = record["id"][len("gene:"):]

        if organism == "human":
            name = record["hs_gene_symbol"]
            down_genes = _get_genes(record, "HGNC", "down_genes")
            up_genes = _get_genes(record, "HGNC", "up_genes")
        elif organism == "rat":
            name = record.get("rs_gene_symbol")
            if name is None:
                # they did not include consistent curation of rat gene symbols
                name = _rat_name_from_human_name(record["hs_gene_symbol"])
            if name is None:
                return
            down_genes = _get_genes(record, "RGD", "down_genes")
            up_genes = _get_genes(record, "RGD", "up_genes")
        elif organism == "mouse":
            name = record["mm_gene_symbol"]
            down_genes = _get_genes(record, "MGI", "down_genes")
            up_genes = _get_genes(record, "MGI", "up_genes")
        else:
            raise ValueError(f"unhandled organism: {organism}")

        pert_type = record["pert_type"]
        up_stmt_cls = UP_MAP.get(pert_type)
        down_stmt_cls = DOWN_MAP.get(pert_type)
        if up_stmt_cls is None or down_stmt_cls is None:
            # The perturbation wasn't readily mappable to an INDRA-like statement class
            return

        # TODO how to use the following metadata?
        geo_id = record["geo_id"]
        cell_type = record["cell_type"]
        evidence = Evidence(
            source_api="creeds",
            annotations={
                "organism": organism,
                "cell": cell_type,
            },
        )
        subject = get_standard_agent(name, {"EGID": perturbed_ncbigene_id})
        for prefix, identifier, name in up_genes:
            target = get_standard_agent(name, {prefix: identifier})
            yield up_stmt_cls(subject, target, copy(evidence))
        for prefix, identifier, name in down_genes:
            target = get_standard_agent(name, {prefix: identifier})
            yield down_stmt_cls(subject, target, copy(evidence))


if __name__ == '__main__':
    GeneProcessor.cli()
