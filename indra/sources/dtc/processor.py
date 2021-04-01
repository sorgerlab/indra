# -*- coding: utf-8 -*-

"""Processor for `Drug Target Commons <https://drugtargetcommons.fimm.fi/>`_."""

import csv
import logging
import pystow
from indra.ontology.standardize import get_standard_agent
from indra.statements import Activation, Agent, Evidence, Inhibition, Statement
from typing import Iterable, List, Optional, Type

logger = logging.getLogger(__name__)

DTC_URL = 'https://drugtargetcommons.fimm.fi/static/Excell_files/DTC_data.csv'

INCREASE_ACTIVITY_RELATIONS = {
    'activation',
}
DECREASE_ACTIVITY_RELATIONS = {
    'inhibition',
    'growth_inhibition',
    'inverse_agonist',
    'cytotoxocity',  # Checked manually and they are inhibitors
}

_UNHANDLED = set()


def _get_statement_cls(s) -> Optional[Type[Statement]]:
    if s in INCREASE_ACTIVITY_RELATIONS:
        return Activation
    if s in DECREASE_ACTIVITY_RELATIONS:
        return Inhibition
    if s in _UNHANDLED:
        _UNHANDLED.add(s)
        logger.warning('unhandled relation: %s', s)


# TODO should names be required for chemicals?
NECESSARY_COLUMNS = [
    # 0,  # compound_id
    # 4,  # target id
    15,  # ep_action_mode
]


def get_rows():
    path = pystow.ensure('indra', 'sources', 'dtc', url=DTC_URL)
    with open(path) as file:
        reader = csv.reader(file, quotechar='"', quoting=csv.QUOTE_MINIMAL, delimiter=',')
        header = next(reader)
        for row in reader:
            yv = dict(zip(header, row))
            if any(not row[i] for i in NECESSARY_COLUMNS):
                continue
            yield yv


def process_row(row) -> Optional[Statement]:
    statement_cls = _get_statement_cls(row['ep_action_mode'])
    target = _get_target(row)
    chemical = _get_chemical(row)
    evidence = _get_evidence(row)
    if statement_cls is None or target is None or chemical is None:
        return None
    return statement_cls(chemical, target, evidence=evidence)


def _get_evidence(row) -> List[Evidence]:
    pubmed_id = row['pubmed_id'] or None
    annotations = {}
    for k in ['assaytype', 'assay_cell_line']:
        v = row[k]
        if v:
            annotations[k] = v
    return [Evidence(source_api='dtc', pmid=pubmed_id, annotations=annotations)]


def _get_target(row) -> Optional[Agent]:
    db_refs = {}
    uniprot_id = row['target_id']
    if uniprot_id:
        db_refs['UP'] = uniprot_id
    else:
        return  # TODO there is still a name or list of names that could be salvaged
    if ',' in uniprot_id:
        logger.warning('unhandled list of targets: %s', unipro)
        return
    return get_standard_agent(
        name=row['target_pref_name'],
        db_refs=db_refs,
    )


def _get_chemical(row) -> Optional[Agent]:
    db_refs = {}
    chembl_id = row['compound_id']
    if chembl_id:
        db_refs['CHEMBL'] = chembl_id
    inchikey = row['standard_inchi_key']
    if inchikey:
        db_refs['INCHIKEY'] = inchikey
    name = row['compound_name']
    if not name:
        return
    return get_standard_agent(
        name=name,
        db_refs=db_refs,
    )


def yield_statements() -> Iterable[Statement]:
    """Clean and map database dump."""
    for row in get_rows():
        statement = process_row(row)
        if statement:
            yield statement


def main():
    from tqdm import tqdm

    statements = yield_statements()
    statements = tqdm(statements, desc='DTC', unit_scale=True)
    for _, statement in zip(range(100), statements):
        print(statement)


if __name__ == '__main__':
    main()
