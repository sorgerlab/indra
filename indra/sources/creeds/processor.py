# -*- coding: utf-8 -*-

"""Processors for CREEDS data."""

from copy import copy
from typing import Any, ClassVar, Iterable, List, Mapping, Optional, Tuple, Type

import click
import pandas as pd
import pystow
from tqdm import tqdm

from indra import statements
from indra.databases import hgnc_client
from indra.ontology.bio import bio_ontology
from indra.ontology.standardize import get_standard_agent
from indra.sources.utils import Processor
from indra.statements import Agent, Evidence, Statement

__all__ = [
    "CREEDSGeneProcessor",
    "CREEDSChemicalProcessor",
    "CREEDSDiseaseProcessor",
]

CREEDS_MODULE = pystow.module("bio", "creeds")
BASE_URL = "http://amp.pharm.mssm.edu/CREEDS/download"
GENE_DATA_URL = f"{BASE_URL}/single_gene_perturbations-v1.0.json"
DISEASE_DATA_URL = f"{BASE_URL}/disease_signatures-v1.0.json"
CHEMICAL_DATA_URL = f"{BASE_URL}/single_drug_perturbations-v1.0.json"

#: A mapping from labels used in CREEDS for species to their
#: organism-specific nomenclature name
ORGANISMS_TO_NS = {
    "mouse": "MGI",
    "human": "HGNC",
    "rat": "RGD",
}
#: A mapping of strings used in the "pert_type" entries in
#: CREEDS data to normalized keys. Several are not curated
#: because they do not readily map to an INDRA statement
PERTURBATIONS = {
    # knockout
    "ko": "knockout",
    "deletion": "knockout",
    "null mutation": "knockout",
    # knockdown
    "kd": "knockdown",
    "deficiency (mutation)": "knockdown",
    "silencing": "knockdown",
    "heterozygotic knockout (nf1+/-)": "knockout",
    # increase
    "induction": "increase",
    "knock-in": "increase",
    "oe": "increase",
    "overexpression": "increase",
    "stimulation of gene product": "increase",
    # activation
    "agonist activation": "activation",
    "drugactivation": "activation",
    "activemutant": "activation",
    "activation (deltanb-cateniner transgenics)": "activation",
    # inhibition
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
    "knockout": statements.DecreaseAmount,
    "knockdown": statements.DecreaseAmount,
    "inhibition": statements.DecreaseAmount,
    "increase": statements.IncreaseAmount,
    "activation": statements.IncreaseAmount,
}
#: A mapping of perturbation types to statement types for the target
#: genes whose gene expression decreases from the given perturbation
#: of the subject gene
DOWN_MAP: Mapping[str, Type[statements.RegulateAmount]] = {
    "knockout": statements.IncreaseAmount,
    "knockdown": statements.IncreaseAmount,
    "inhibition": statements.IncreaseAmount,
    "increase": statements.DecreaseAmount,
    "activation": statements.DecreaseAmount,
}


def _process_pert_type(s: str) -> str:
    x = s.strip().lower()
    return PERTURBATIONS.get(x, x)



def _get_genes(
    record: Mapping[str, Any],
    prefix: str,
    key: str,
) -> List[Tuple[str, str, str]]:
    rv = []
    #: A list of 2-tuples with the gene symbol then the expression value
    expressions = record[key]
    for symbol, _ in expressions:
        try:
            # FIXME this doesn't currently succeed for MGI nor RGD
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


def _get_evidence(record: Mapping[str, Any]) -> Evidence:
    # TODO how to use the following metadata?
    geo_id = record["geo_id"]
    cell_type = record["cell_type"]
    organism = record["organism"]
    return Evidence(
        source_api="creeds",
        annotations={
            "organism": organism,
            "cell": cell_type,
            "geo": geo_id,
        },
    )


def _get_regulations(
    record: Mapping[str, Any],
) -> Tuple[List[Tuple[str, str, str]], List[Tuple[str, str, str]]]:
    organism = record["organism"]
    prefix = ORGANISMS_TO_NS[organism]
    up_genes = _get_genes(record, prefix, "up_genes")
    down_genes = _get_genes(record, prefix, "down_genes")
    return up_genes, down_genes


def _process_record_helper(
    record, subject, up_stmt_cls, down_stmt_cls
) -> Iterable[Statement]:
    up_genes, down_genes = _get_regulations(record)
    evidence = _get_evidence(record)
    for prefix, identifier, name in up_genes:
        target = get_standard_agent(name, {prefix: identifier})
        yield up_stmt_cls(subject, target, copy(evidence))
    for prefix, identifier, name in down_genes:
        target = get_standard_agent(name, {prefix: identifier})
        yield down_stmt_cls(subject, target, copy(evidence))


class CREEDSProcessor(Processor):
    """A base processor for CREEDS, which can either take records directly,
    or looks for it on the web using the given class variable ``url`` in
    combination with :mod:`pystow`.
    """

    #: The URL of the remote JSON file
    url: ClassVar[str]
    #: The processed statements (after ``extract_statements()`` is run)
    statements: List[Statement]

    def __init__(self, records: Optional[List[Mapping[str, Any]]] = None):
        self.records = records or CREEDS_MODULE.ensure_json(url=self.url)
        self.statements = []

    def extract_statements(self) -> List[Statement]:
        """Generate and store statements if not pre-cached, then return then."""
        if not self.statements:
            self.statements = list(self.iter_statements())
        return self.statements

    def iter_statements(self) -> Iterable[Statement]:
        for record in tqdm(self.records, desc=f"Processing {self.name}"):
            yield from self.process_record(record)

    @classmethod
    def process_record(cls, record: Mapping[str, Any]) -> Iterable[Statement]:
        raise NotImplementedError


LOGGED_MISSING_PART = set()


class CREEDSGeneProcessor(CREEDSProcessor):
    """A processor for single gene perturbation experiments in CREEDS."""

    name = "creeds_gene"
    url = GENE_DATA_URL

    @staticmethod
    def get_subject(record) -> Optional[Agent]:
        perturbed_ncbigene_id = record["id"][len("gene:") :]
        organism = record["organism"]
        if organism == "human":
            name = record["hs_gene_symbol"]
        elif organism == "rat":
            name = record.get("rs_gene_symbol")
            if name is None:
                # they did not include consistent curation of rat gene symbols
                name = _rat_name_from_human_name(record["hs_gene_symbol"])
            if name is None:
                return
        elif organism == "mouse":
            name = record["mm_gene_symbol"]
        else:
            raise ValueError(f"unhandled organism: {organism}")
        return get_standard_agent(name, {"EGID": perturbed_ncbigene_id})

    @classmethod
    def process_record(cls, record: Mapping[str, Any]) -> Iterable[Statement]:
        subject = cls.get_subject(record)
        if subject is None:
            return

        pert_type = _process_pert_type(record["pert_type"])
        up_stmt_cls = UP_MAP.get(pert_type)
        down_stmt_cls = DOWN_MAP.get(pert_type)
        if up_stmt_cls is None or down_stmt_cls is None:
            if pert_type not in LOGGED_MISSING_PART:
                tqdm.write(f"Could not look up pert_type {pert_type}")
                LOGGED_MISSING_PART.add(pert_type)
            return

        yield from _process_record_helper(record, subject, up_stmt_cls, down_stmt_cls)


class CREEDSDiseaseProcessor(CREEDSProcessor):
    """A processor for disease perturbation experiments in CREEDS."""

    name = "creeds_disease"
    url = DISEASE_DATA_URL

    @staticmethod
    def get_subject(record) -> Agent:
        xrefs = {}
        doid = record["do_id"]
        if doid:
            xrefs["DOID"] = doid
        umls_id = record["umls_cui"]
        if umls_id:
            xrefs["UMLS"] = umls_id
        name = record["disease_name"]
        return get_standard_agent(name, xrefs)

    @classmethod
    def process_record(cls, record) -> Iterable[Statement]:
        subject = cls.get_subject(record)
        yield from _process_record_helper(
            record,
            subject,
            up_stmt_cls=statements.IncreaseAmount,
            down_stmt_cls=statements.DecreaseAmount,
        )


class CREEDSChemicalProcessor(CREEDSProcessor):
    """A processor for chemical perturbation experiments in CREEDS."""

    name = "creeds_chemical"
    url = CHEMICAL_DATA_URL

    @staticmethod
    def get_subject(record) -> Agent:
        xrefs = {}
        smiles = record["smiles"]
        if smiles:
            xrefs["SMILES"] = smiles
        pubchem_compound_id = record["pubchem_cid"]
        if pubchem_compound_id:
            xrefs["PUBCHEM"] = str(pubchem_compound_id)
        drugbank_id = record["drugbank_id"]
        if drugbank_id:
            xrefs["DRUGBANK"] = drugbank_id
        name = record["drug_name"]
        return get_standard_agent(name, xrefs)

    @classmethod
    def process_record(cls, record) -> Iterable[Statement]:
        subject = cls.get_subject(record)
        yield from _process_record_helper(
            record,
            subject,
            up_stmt_cls=statements.IncreaseAmount,
            down_stmt_cls=statements.DecreaseAmount,
        )


@click.command()
@click.pass_context
def _main(ctx: click.Context):
    ctx.invoke(CREEDSChemicalProcessor.get_cli())
    ctx.invoke(CREEDSGeneProcessor.get_cli())
    ctx.invoke(CREEDSDiseaseProcessor.get_cli())


if __name__ == "__main__":
    _main()
