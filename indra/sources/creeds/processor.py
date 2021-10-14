from collections import Counter
from typing import Iterable

import pandas as pd
import pystow
from tabulate import tabulate

from indra.ontology.bio import bio_ontology
from indra.ontology.standardize import get_standard_agent
from indra.sources.utils import Processor

__all__ = [
    "GeneProcessor",
]

from indra.statements import Statement

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


def preprocess_metadata(force: bool = False) -> pd.DataFrame:
    metadata = CREEDS_MODULE.ensure_csv(url=GENE_PERTURBATIONS_METADATA_URL, force=force, read_csv_kwargs=dict(sep=","))
    metadata = metadata[metadata.pert_type.notna()]
    metadata.id = metadata.id.map(lambda s: s[len("gene:"):])
    metadata.organism = metadata.id.map(ORGANISMS, na_action="ignore")

    metadata.pert_type = metadata.pert_type.map(str.strip)
    metadata.pert_type = metadata.pert_type.map(str.lower)
    for key in metadata.pert_type.unique():
        if key not in PERTURBATIONS and pd.notna(key):
            PERTURBATIONS[key] = key
    metadata.pert_type = metadata.pert_type.map(PERTURBATIONS)
    return metadata


def preprocess_records() -> Iterable[Statement]:
    records = CREEDS_MODULE.ensure_json(url=GENE_PERTURBATIONS_DATA_URL)
    for record in records:
        yield from process_record(record)


def x(record, ns, key):
    return [
        (ns, symbol, bio_ontology.get_id_from_name(ns, symbol))
        for symbol, _ in record[key]
    ]


def process_record(record):
    organism = record["organism"]
    perturbed_ncbigene_id = record["id"][len("gene:"):]

    if organism == "human":
        name = record["hs_gene_symbol"]
        down_genes = x(record, "HGNC", "down_genes")
        up_genes = x(record, "HGNC", "up_genes")
    elif organism == "rat":
        name = record["rs_gene_symbol"]
        down_genes = x(record, "RGD", "down_genes")
        up_genes = x(record, "RGD", "up_genes")
    elif organism == "mouse":
        name = record["mm_gene_symbol"]
        down_genes = x(record, "MGI", "down_genes")
        up_genes = x(record, "MGI", "up_genes")
    else:
        raise ValueError
    subject = get_standard_agent(name, {"EGID": perturbed_ncbigene_id})
    stmt_cls = ...
    for ns, name, id in up_genes:
        if id is None:
            continue
        target = get_standard_agent(name, {ns: id})
        yield stmt_cls(subject, target)

class GeneProcessor(Processor):
    def __init__(self):
        self.metadata = preprocess_metadata()

        print(tabulate(Counter(self.metadata.pert_type).most_common(25)))

        self.records = {}

    def extract_statements(self):
        return []

    @staticmethod
    def _process_record(record):
        pass


if __name__ == '__main__':
    GeneProcessor.cli()
