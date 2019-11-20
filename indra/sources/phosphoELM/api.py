import csv

ppelm_s3_key = ''


def process_from_dump(fname=None):
    ppelm_json = []
    if fname is None:
        # ToDo Get from S3
        pass
    else:
        with open(fname, 'r') as f:
            csv_reader = csv.reader(f.readlines(), dialect='csv')
            columns = next(csv_reader)
            for entry in csv_reader:
                row_dict = {columns[n]: entry[n]
                            for n in range(len(columns))}
                ppelm_json.append(row_dict)
    return ppelm_json
