import io
import pandas
import requests

from .processor import TrrustProcessor


trrust_human_url = 'https://www.grnpedia.org/trrust/data/trrust_rawdata' \
                   '.human.tsv'


def process_from_web():
    res = requests.get(trrust_human_url)
    table = res.read()
    df = pandas.read_table(io.StringIO(table))
