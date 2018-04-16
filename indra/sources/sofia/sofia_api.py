import openpyxl
from .processor import SofiaProcessor

def process_table(fname, sheet_name):
    """Return processor by processing a given sheet of a spreadsheet file.

    Parameters
    ----------
    fname : str
        The name of the Excel file (typically .xlsx extension) to process
    sheet_name : str
        The name of the sheet in the Excel file that has extractions

    Returns
    -------
    sp : indra.sources.sofia.processor.SofiaProcessor
        A SofiaProcessor object which has a list of extracted INDRA
        Statements as its statements attribute
    """
    book = openpyxl.load_workbook(fname, read_only=True)
    sheet = book[sheet_name]
    rows = sheet.rows
    sp = SofiaProcessor(rows)
    return sp
