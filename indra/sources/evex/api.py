import requests
import tqdm

url = 'http://evexdb.org/download/network-format/Metazoa/Homo_sapiens.tar.gz'
standoff_root = 'http://evexdb.org/download/standoff-annotation/version-0.1/'


def download_evex():
    """Download EVEX standoff output."""
    import pystow
    from bs4 import BeautifulSoup
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
