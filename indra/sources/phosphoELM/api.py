import csv

ppelm_s3_key = ''


def process_from_dump(fname=None, delimiter='\t'):
    if fname is None:
        # ToDo Get from S3
        return []
    else:
        with open(fname, 'r') as f:
            csv_reader = csv.reader(f.readlines(), delimiter=delimiter)
            ppelm_json = _get_json_from_entry_rows(csv_reader)
    return ppelm_json


def _get_json_from_entry_rows(row_iter):
    ppelm_json = []
    columns = next(row_iter)
    for entry in row_iter:
        row_dict = {columns[n]: entry[n]
                    for n in range(len(columns))}
        ppelm_json.append(row_dict)
    return ppelm_json
