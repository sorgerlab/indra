import openpyxl
from .processor import SofiaProcessor


def process_table(fname):
    """Return processor by processing a given sheet of a spreadsheet file.

    Parameters
    ----------
    fname : str
        The name of the Excel file (typically .xlsx extension) to process

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute
    """
    book = openpyxl.load_workbook(fname, read_only=True)
    try:
        rel_sheet = book['Relations']
    except Exception as e:
        rel_sheet = book['Causal']
    event_sheet = book['Events']
    entities_sheet = book['Entities']
    sp = SofiaProcessor(rel_sheet.rows, event_sheet.rows, entities_sheet.rows)
    return sp
