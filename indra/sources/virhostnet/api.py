import pandas
from .processor import VirhostnetProcessor

vhn_url = ('http://virhostnet.prabi.fr:9090/psicquic/webservices/current/'\
           'search/query/*')


data_columns = [
    'host_grounding', 'vir_grounding', 'host_mnemonic', 'vir_mnemonic',
    'host_mnemonic2', 'vir_mnemonic2', 'exp_method',
    'not_sure', 'publication', 'host_tax', 'vir_tax',
    'int_type', 'source', 'source_id', 'score'
]


def process_from_web():
    df = pandas.read_csv(vhn_url, delimiter='\t', names=data_columns,
                         header=None)
    return process_df(df)


def process_tsv(fname):
    df = pandas.read_csv(fname, delimiter='\t', names=data_columns,
                         header=None)
    return process_df(df)


def process_df(df):
    vp = VirhostnetProcessor(df)
    vp.extract_statements()
    return vp
