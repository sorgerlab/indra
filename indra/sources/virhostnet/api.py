__all__ = ['process_from_web', 'process_tsv', 'process_df']

import pandas
import logging
from .processor import VirhostnetProcessor

logger = logging.getLogger(__name__)


vhn_url = ('http://virhostnet.prabi.fr:9090/psicquic/webservices/current/'\
           'search/query/')


data_columns = [
    'host_grounding', 'vir_grounding', 'host_mnemonic', 'vir_mnemonic',
    'host_mnemonic2', 'vir_mnemonic2', 'exp_method',
    'dash', 'publication', 'host_tax', 'vir_tax',
    'int_type', 'source', 'source_id', 'score'
]


def process_from_web(query=None):
    """Process host-virus interactions from the VirHostNet website.

    Parameters
    ----------
    query : Optional[str]
        A query that constrains the results to a given subset of the VirHostNet
        database. Example: "taxid:2697049" to search for interactions for
        SARS-CoV-2. If not provided, By default, the "*" query is used which
        returns the full database.

    Returns
    -------
    VirhostnetProcessor
        A VirhostnetProcessor object which contains a list of extracted
        INDRA Statements in its statements attribute.
    """
    # Search for everything to get the full download by default
    url = vhn_url + ('*' if query is None else query)
    logger.info('Processing VirHostNet data from %s' % url)
    df = pandas.read_csv(url, delimiter='\t', names=data_columns,
                         header=None)
    return process_df(df)


def process_tsv(fname):
    """Process a TSV data file obtained from VirHostNet.

    Parameters
    ----------
    fname : str
        The path to the VirHostNet tabular data file (in the same format as
        the web service).

    Returns
    -------
    VirhostnetProcessor
        A VirhostnetProcessor object which contains a list of extracted
        INDRA Statements in its statements attribute.
    """
    df = pandas.read_csv(fname, delimiter='\t', names=data_columns,
                         header=None)
    return process_df(df)


def process_df(df):
    """Process a VirHostNet pandas DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame representing VirHostNet interactions (in the same format as
        the web service).

    Returns
    -------
    VirhostnetProcessor
        A VirhostnetProcessor object which contains a list of extracted
        INDRA Statements in its statements attribute.
    """
    vp = VirhostnetProcessor(df)
    vp.extract_statements()
    return vp
