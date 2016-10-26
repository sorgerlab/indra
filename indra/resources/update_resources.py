import os
import logging
import shutil
import requests

path = os.path.join(os.path.dirname(__file__), 'tmp')
logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('update_resources')
logger.setLevel(logging.INFO)

def update_hgnc_entries():
    url = 'http://tinyurl.com/gnv32vh'
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    fname = os.path.join(path, 'hgnc_entries.tsv')
    logger.info('Writing into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(res.content)

if __name__ == '__main__':
    update_hgnc_entries()
