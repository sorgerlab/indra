import os
import glob
import logging
import pickle
import tarfile
import requests
import pandas
import pystow
import tqdm

from .processor import EvexProcessor

logger = logging.getLogger(__name__)

human_network = 'http://evexdb.org/download/network-format/Metazoa/' \
    'Homo_sapiens.tar.gz'
standoff_root = 'http://evexdb.org/download/standoff-annotation/version-0.1/'


def process_human_events():
    standoff_index = build_standoff_index()
    network_file = pystow.ensure('evex', name='Homo_sapiens.tar.gz',
                                 url=human_network)
    with tarfile.open(network_file, 'r:gz') as fh:
        relations_file = fh.extractfile('EVEX_relations_9606.tab')
        articles_file = fh.extractfile('EVEX_articles_9606.tab')
        relations_df = pandas.read_csv(relations_file, sep='\t')
        articles_df = pandas.read_csv(articles_file, sep='\t')
    ep = EvexProcessor(relations_df, articles_df, standoff_index)
    ep.process_statements()
    return ep


def build_standoff_index(cached=True):
    """Build an index of publications in standoff bulk archive files."""
    cache_file = pystow.join('evex', name='standoff_index.pkl')
    if cached and cache_file.exists():
        logger.info('Loading standoff index from %s' % cache_file.as_posix())
        with open(cache_file, 'rb') as fh:
            return pickle.load(fh)
    index = {}
    for fname in tqdm.tqdm(glob.glob(os.path.join(
                                     pystow.join('evex').as_posix(), 'batch*')),
                           desc='Building standoff index'):
        try:
            with tarfile.open(fname, 'r:gz') as fh:
                names = fh.getnames()
        except tarfile.ReadError:
            logger.error('Could not read tarfile %s' % fname)
            continue
        ids = {tuple(os.path.splitext(name)[0].split('_')[:2])
               for name in names if name.endswith('ann')}
        for paper_id in ids:
            index[paper_id] = fname
    if cached:
        with open(cache_file, 'wb') as fh:
            pickle.dump(index, fh)
    return index


def download_evex():
    """Download EVEX standoff output."""
    from bs4 import BeautifulSoup
    # Download human network first
    pystow.ensure('evex', name='Homo_sapiens.tar.gz', url=human_network)
    # Now download all the standoff files
    res = requests.get(standoff_root)
    soup = BeautifulSoup(res.text, 'html.parser')
    children = [standoff_root + node.get('href')
                for node in soup.find_all('a')
                if node.get('href').startswith('files')]
    for child in tqdm.tqdm(children):
        res = requests.get(child)
        soup = BeautifulSoup(res.text, 'html.parser')
        downloadables = [child + node.get('href')
                         for node in soup.find_all('a')
                         if node.get('href').startswith('batch')]
        for downloadable in downloadables:
            fname = downloadable.split('/')[-1]
            pystow.ensure('evex', name=fname, url=downloadable)
