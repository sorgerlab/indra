import openpyxl
from .processor import SofiaProcessor

def process_table(fname, sheet_name):
    book = openpyxl.load_workbook(fname, read_only=True)
    sheet = book[sheet_name]
    rows = sheet.rows
    sp = SofiaProcessor(rows)
    return sp
