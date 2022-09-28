import os
import glob
import logging
import pickle
import tarfile
from urllib.request import urlretrieve
import requests
import pandas
import tqdm

from .processor import EvexProcessor

logger = logging.getLogger(__name__)

human_network = 'http://evexdb.org/download/network-format/Metazoa/' \
    'Homo_sapiens.tar.gz'
standoff_root = 'http://evexdb.org/download/standoff-annotation/version-0.1/'


def process_human_events(base_folder=None):
    """Process all human events available in EVEX.

    Note that unless the standoff files have already been downloaded using the
    `download_evex` function, the Statements produced by this function
    will not carry any evidence text, agent text and various other metadata
    in them for which the standoff files are required.

    Parameters
    ----------
    base_folder : Optional[str]
        If provided, the given base folder is used to download the human
        network file from EVEX. Otherwise, the `pystow` package is used
        to create an `evex` folder within the pystow base path,
        typically ~/.data/evex.

    Returns
    -------
    EvexProcessor
        An EvexProcessor instance with the extracted INDRA Statements
        as its statements attribute.
    """
    if not base_folder:
        import pystow
        base_folder = pystow.join('evex').as_posix()
    standoff_index = build_standoff_index()
    network_file = os.path.join(base_folder, 'Homo_sapiens.tar.gz')
    if not os.path.exists(network_file):
        urlretrieve(human_network, network_file)
    with tarfile.open(network_file, 'r:gz') as fh:
        relations_file = fh.extractfile('EVEX_relations_9606.tab')
        articles_file = fh.extractfile('EVEX_articles_9606.tab')
        relations_df = pandas.read_csv(relations_file, sep='\t')
        articles_df = pandas.read_csv(articles_file, sep='\t')
    ep = EvexProcessor(relations_df, articles_df, standoff_index)
    ep.process_statements()
    return ep


def build_standoff_index(cached=True, base_folder=None):
    """Build an index of publications in standoff bulk archive files.

    This index is necessary to figure out which standoff archive the annotations
    for a given article are in.

    Parameters
    ----------
    cached: Optional[bool]
        If True, the standoff index is cached in the base folder and isn't
        regenerated if this function is called again, just reloaded.
        This is useful since generating the full standoff file index
        can take a long time. Default: True
    base_folder : Optional[str]
        If provided, the given base folder is used to download the human
        network file from EVEX. Otherwise, the `pystow` package is used
        to create an `evex` folder within the pystow base path,
        typically ~/.data/evex.
    """
    if not base_folder:
        import pystow
        base_folder = pystow.join('evex').as_posix()
    cache_file = os.path.join(base_folder, 'standoff_index.pkl')
    if cached and os.path.exists(cache_file):
        logger.info('Loading standoff index from %s' % cache_file)
        with open(cache_file, 'rb') as fh:
            return pickle.load(fh)
    index = {}
    for fname in tqdm.tqdm(glob.glob(os.path.join(base_folder, 'batch*')),
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


def download_evex(base_folder=None):
    """Download EVEX human network and standoff output files.

    This function downloads the human network file as well as a large number
    of standoff output files. These files are necessary to find evidence text,
    agent text and agent coordinates to be used in INDRA. Note that there
    are over 4 thousand such files, and the overall size is around 6 GB.

    Parameters
    ----------
    base_folder : Optional[str]
        If provided, the given base folder is used to download the human
        network file from EVEX. Otherwise, the `pystow` package is used
        to create an `evex` folder within the pystow base path,
        typically ~/.data/evex.
    """
    from bs4 import BeautifulSoup
    if not base_folder:
        import pystow
        base_folder = pystow.join('evex').as_posix()
    # Download human network first
    fname = os.path.join(base_folder, 'Homo_sapiens.tar.gz')
    if not os.path.exists(fname):
        urlretrieve(human_network, fname)
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
            fname = os.path.join(base_folder, downloadable.split('/')[-1])
            if not os.path.exists(fname):
                urlretrieve(downloadable, fname)
