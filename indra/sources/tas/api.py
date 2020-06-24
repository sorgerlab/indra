from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

__all__ = ['process_csv', 'process_from_web']

import csv
import logging
import requests
from hashlib import md5

from .processor import TasProcessor
from indra.util import read_unicode_csv

tas_data_url = 'https://bigmech.s3.amazonaws.com/indra-db/tas.csv'
tas_resource_md5 = '554ccba4617aae7b3b06a62893424c7f'


logger = logging.getLogger(__name__)


def _load_data(data_iter):
    # Get the headers.
    headers = data_iter[0]

    # For some reason this heading is oddly formatted and inconsistent with the
    # rest, or with the usual key-style for dicts.
    data = [{header: val for header, val in zip(headers, line)}
            for line in data_iter[1:]]
    return data


def process_from_web(affinity_class_limit=2):
    """Return a TasProcessor for the contents of the TAS dump online.

    Interactions are classified into the following classes based on affinity:
      | 1  -- Kd < 100nM
      | 2  -- 100nM < Kd < 1uM
      | 3  -- 1uM < Kd < 10uM
      | 10 -- Kd > 10uM
    By default, only classes 1 and 2 are extracted but the affinity_class_limit
    parameter can be used to change the upper limit of extracted classes.

    Parameters
    ----------
    affinity_class_limit : Optional[int]
        Defines the highest class of binding affinity that is included in the
        extractions. Default: 2

    Returns
    -------
    TasProcessor
        A TasProcessor object which has a list of INDRA Statements extracted
        from the CSV file representing drug-target inhibitions in its
        statements attribute.
    """
    logger.info('Downloading TAS data from %s' % tas_data_url)
    res = requests.get(tas_data_url)
    observed_checksum = md5(res.text.encode('utf-8')).hexdigest()
    logger.info('Verifying md5 checksum of data')
    if tas_resource_md5 != observed_checksum:
        raise RuntimeError('Checksum for downloaded TAS data does not'
                           ' match expected value')
    res.raise_for_status()
    logger.info('Finished downloading TAS data from %s' % tas_data_url)
    data_iter = list(csv.reader(res.text.splitlines(), delimiter=','))
    return TasProcessor(_load_data(data_iter), affinity_class_limit)


def process_csv(fname, affinity_class_limit=2):
    """Return a TasProcessor for the contents of a given CSV file..

    Interactions are classified into the following classes based on affinity:
      | 1  -- Kd < 100nM
      | 2  -- 100nM < Kd < 1uM
      | 3  -- 1uM < Kd < 10uM
      | 10 -- Kd > 10uM
    By default, only classes 1 and 2 are extracted but the affinity_class_limit
    parameter can be used to change the upper limit of extracted classes.

    Parameters
    ----------
    fname : str
        The path to a local CSV file containing the TAS data.
    affinity_class_limit : Optional[int]
        Defines the highest class of binding affinity that is included in the
        extractions. Default: 2

    Returns
    -------
    TasProcessor
        A TasProcessor object which has a list of INDRA Statements extracted
        from the CSV file representing drug-target inhibitions in its
        statements attribute.
    """
    data_iter = list(read_unicode_csv(fname))
    return TasProcessor(_load_data(data_iter), affinity_class_limit)
