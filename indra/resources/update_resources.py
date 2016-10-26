import os
import shutil
import urllib
import urllib2
import logging
import requests
from contextlib import closing

path = os.path.join(os.path.dirname(__file__), 'tmp')
logging.basicConfig(format='%(levelname)s: indra/%(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('update_resources')
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('requests').setLevel(logging.ERROR)
logger.setLevel(logging.INFO)

def update_hgnc_entries():
    logger.info('----------------------------')
    logger.info('--Updating HGNC entries-----')
    url = 'http://tinyurl.com/gnv32vh'
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    fname = os.path.join(path, 'hgnc_entries.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(res.content)
    logger.info('----------------------------')

def update_kinases():
    logger.info('----------------------------')
    logger.info('--Updating kinase list------')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=entry_name&desc=no&compress=no&query=database:(type:' + \
        'interpro%20ipr000719)%20AND%20reviewed:yes%20AND%20organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&fil=&force=no' + \
        '&format=tab&columns=id,genes(PREFERRED),organism-id,entry%20name'
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    fname = os.path.join(path, 'kinases.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(res.content)
    logger.info('----------------------------')

def update_uniprot_entries():
    logger.info('----------------------------')

    logger.info('--Updating UniProt entries--')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'fil=&force=no&format=tab&columns=id,genes(PREFERRED),organism-id,' + \
        'entry%20name'
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download %s' % url)
        return
    fname = os.path.join(path, 'uniprot_entries.tsv')
    logger.info('Saving into %s' % fname)
    with open(fname, 'wb') as fh:
        fh.write(res.content)
    logger.info('----------------------------')

def update_uniprot_sec_ac():
    logger.info('----------------------------')
    logger.info('--Updating UniProt secondary accession--')
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
        'docs/sec_ac.txt'
    logger.info('Downloading %s' % url)
    fname = os.path.join(path, 'uniprot_sec_ac.txt')
    with closing(urllib.urlopen(url)) as res:
        with open(fname, 'wb') as fh:
            logger.info('Saving into %s' % fname)
            shutil.copyfileobj(res, fh)
    with open(fname, 'wb') as fh:
        fh.write(res.content)
    logger.info('----------------------------')

def update_uniprot_subcell_loc():
    logger.info('----------------------------')
    logger.info('--Updating UniProt secondary accession--')
    url = 'http://www.uniprot.org/locations/?' + \
        '%20sort=&desc=&compress=no&query=&force=no&format=tab&columns=id'
    logger.info('Downloading %s' % url)
    fname = os.path.join(path, 'uniprot_subcell_loc.tsv')
    with closing(urllib.urlopen(url)) as res:
        with open(fname, 'wb') as fh:
            logger.info('Saving into %s' % fname)
            shutil.copyfileobj(res, fh)
    with open(fname, 'wb') as fh:
        fh.write(res.content)
    logger.info('----------------------------')

if __name__ == '__main__':
    #update_hgnc_entries()
    #update_kinases()
    #update_uniprot_entries()
    #update_uniprot_sec_ac()
    update_uniprot_subcell_loc()
