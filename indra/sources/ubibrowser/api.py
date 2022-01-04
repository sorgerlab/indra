import pandas
from .processor import UbiBrowserProcessor


DOWNLOAD_URL = 'http://ubibrowser.ncpsb.org.cn/v2/Public/download/literature/'
E3_URL = DOWNLOAD_URL + 'literature.E3.txt'
DUB_URL = DOWNLOAD_URL + 'literature.DUB.txt'


def process_from_web():
    e3_df = pandas.read_csv(E3_URL, sep='\t')
    dub_df = pandas.read_csv(DUB_URL, sep='\t')
    return process_df(e3_df, dub_df)


def process_df(e3_df, dub_df):
    up = UbiBrowserProcessor(e3_df, dub_df)
    up.extract_statements()
    return up