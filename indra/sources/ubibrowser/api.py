__all__ = ['process_from_web', 'process_df']

import pandas
from .processor import UbiBrowserProcessor


DOWNLOAD_URL = 'http://ubibrowser.ncpsb.org.cn/v2/Public/download/literature/'
E3_URL = DOWNLOAD_URL + 'literature.E3.txt'
DUB_URL = DOWNLOAD_URL + 'literature.DUB.txt'


def process_from_web() -> UbiBrowserProcessor:
    """Download the UbiBrowser data from the web and process it.

    Returns
    -------
    :
        An UbiBrowserProcessor object with INDRA Statements
        extracted in its statements attribute.
    """
    e3_df = pandas.read_csv(E3_URL, sep='\t')
    dub_df = pandas.read_csv(DUB_URL, sep='\t')
    return process_df(e3_df, dub_df)


def process_df(e3_df: pandas.DataFrame, dub_df: pandas.DataFrame) \
        -> UbiBrowserProcessor:
    """Process data frames containing UbiBrowser data.

    Parameters
    ----------
    e3_df :
        A data frame containing UbiBrowser E3 data.
    dub_df :
        A data frame containing UbiBrowser DUB data.

    Returns
    -------
    :
        An UbiBrowserProcessor object with INDRA Statements
        extracted in its statements attribute.
    """
    up = UbiBrowserProcessor(e3_df, dub_df)
    up.extract_statements()
    return up