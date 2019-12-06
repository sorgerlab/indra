import csv
import logging

from .processor import PhosphoElmProcessor

logger = logging.getLogger(__name__)


def process_from_dump(fname, delimiter='\t'):
    """Process a Phospho.ELM file dump

    The dump can be obtained at http://phospho.elm.eu.org/dataset.html.

    Parameters
    ----------
    fname : str
        File path to the phospho.ELM file dump.
    delimiter : str
        The delimiter to use for csv.reader

    Returns
    -------
    indra.sources.phosphoelm.PhosphoElmProcessor
        An instance of a PhosphoElmProcessor containing the statements
        generated from the file dump
    """
    with open(fname, 'r') as f:
        csv_reader = csv.reader(f, delimiter=delimiter)
        ppelm_json = _get_json_from_entry_rows(csv_reader)
    pep = PhosphoElmProcessor(ppelm_json)
    pep.process_phosphorylations()
    return pep


def _get_json_from_entry_rows(row_iter):
    """Loop body to generate a json friendly structure"""
    ppelm_json = []
    columns = next(row_iter)
    for entry in row_iter:
        row_dict = {c: e for c, e in zip(columns, entry)}
        ppelm_json.append(row_dict)
    return ppelm_json
