import csv

from .processor import PhosphoELMPRocessor

s3_bucket = 'bigmech'
ppelm_s3_key = ''


def process_from_dump(fname=None, delimiter='\t'):
    """Process a phospho.ELM file dump

    fname : str
        File path to the phospho.ELM file dump
    delimiter : str
        The delimiter to use for csv.reader

    Returns
    -------
    ppp : indra.sources.phosphoELM.PhosphoELMPRocessor
        An instance of a PhosphoELMPRocessor containing the statements
        generated from the file dump
    """
    if fname is None:
        # ToDo Get from S3
        return []
    else:
        with open(fname, 'r') as f:
            csv_reader = csv.reader(f.readlines(), delimiter=delimiter)
            ppelm_json = _get_json_from_entry_rows(csv_reader)
    return PhosphoELMPRocessor(file_dump_json=ppelm_json)


def _get_json_from_entry_rows(row_iter):
    """Loop body to generate a json friendly structure"""
    ppelm_json = []
    columns = next(row_iter)
    for entry in row_iter:
        row_dict = {columns[n]: entry[n]
                    for n in range(len(columns))}
        ppelm_json.append(row_dict)
    return ppelm_json
