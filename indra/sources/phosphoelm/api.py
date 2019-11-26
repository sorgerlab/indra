import csv
import logging

from .processor import PhosphoElmProcessor

logger = logging.getLogger(__name__)

s3_bucket = 'bigmech'
ppelm_s3_key = 'indra-db/external_databases/phosphoELM_all_2015-04.dump'
kinases_list_web = 'http://phospho.elm.eu.org/kinases.html'


def process_from_dump(fname=None, delimiter='\t'):
    """Process a phospho.ELM file dump

    Parameters
    ----------
    fname : str
        File path to the phospho.ELM file dump. If none is provided,
        the file will be downloaded from S3.
    delimiter : str
        The delimiter to use for csv.reader

    Returns
    -------
    pep : indra.sources.phosphoelm.PhosphoElmProcessor
        An instance of a PhosphoElmProcessor containing the statements
        generated from the file dump
    """
    if fname is None:
        s3 = _get_s3_client()
        s3_obj = s3.get_object(Bucket=s3_bucket, Key=ppelm_s3_key)
        csv_reader = csv.reader(
            s3_obj['Body'].read().decode('utf8').splitlines(True),
            delimiter='\t'
        )
    else:
        with open(fname, 'r') as f:
            csv_reader = csv.reader(f.readlines(), delimiter=delimiter)
    ppelm_json = _get_json_from_entry_rows(csv_reader)
    pep = PhosphoElmProcessor(ppelm_json)
    pep.process_phosphorylations()
    return pep


def _get_json_from_entry_rows(row_iter):
    """Loop body to generate a json friendly structure"""
    ppelm_json = []
    columns = next(row_iter)
    for entry in row_iter:
        row_dict = {columns[n]: entry[n]
                    for n in range(len(columns))}
        ppelm_json.append(row_dict)
    return ppelm_json


def _get_s3_client(unsigned=True):
    import boto3
    from botocore import UNSIGNED
    from botocore.client import Config

    """Return a boto3 S3 client with optional unsigned config.

    Parameters
    ----------
    unsigned : Optional[bool]
        If True, the client will be using unsigned mode in which public
        resources can be accessed without credentials. Default: True

    Returns
    -------
    botocore.client.S3
        A client object to AWS S3.
    """
    if unsigned:
        return boto3.client('s3', config=Config(signature_version=UNSIGNED))
    else:
        return boto3.client('s3')
