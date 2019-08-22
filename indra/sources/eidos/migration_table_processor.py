import openpyxl


def process_sheet(sheet):
    for row in sheet.rows:
        for cell in row:
            print(cell)
            import ipdb; ipdb.set_trace()
            break


if __name__ == '__main__':
    fname = 'Initial annotation exercise for migration use case.xlsx'

    wb = openpyxl.load_workbook(fname, read_only=True)
    sheets = wb.sheetnames
    cag_sheets = [s for s in sheets if 'CAG' in s]

    for sheet_name in cag_sheets:
        sheet = wb[sheet_name]
        process_sheet(sheet)
