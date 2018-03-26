from indra.sources.biogrid.processor import BiogridProcessor
import requests
import logging
import tempfile
import zipfile
import os
import shutil

logger = logging.getLogger('cwms_rdf')

biogrid_file_url = 'https://downloads.thebiogrid.org/Download/BioGRID/' + \
        'Release-Archive/BIOGRID-3.4.158/BIOGRID-ALL-3.4.158.tab2.zip'
biogrid_text_file = 'BIOGRID-ALL-3.4.158.tab2.txt'

def process_file(filename=None):
    """Processes a biogrid tab-separated filer.

    Parameters:
    -----------
    text: str
        The filename of the biogrid file, if None downloaded from the web.

    Returns
    -------
    bp: indra.sources.cwms.BiogridProcessor
        A BiogridProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    tmp_dir = None
    if filename == None:
        # Filename not specified, download from the web
        tmp_dir = tempfile.mkdtemp('indra_biogrid')

        target_zip = os.path.join(tmp_dir, 'biogrid.zip')
        _download_file(biogrid_file_url, target_zip)

        zip_file = zipfile.ZipFile(target_zip)
        zip_file.extract(biogrid_text_file, tmp_dir)

        filename = os.path.join(tmp_dir, biogrid_text_file)

    bp = BiogridProcessor(filename)

    if tmp_dir is not None:
        # Clean up temporary files
        shutil.rmtree(tmp_dir)
    return bp

def _download_file(url, target):
    """Downloads a file to the specified target, works in python 2 and
    python 3."""
    r = requests.get(url)
    with open(target, 'wb') as f:
        f.write(r.content)
    
if __name__ == '__main__':
    fname = '/Users/daniel/Downloads/BIOGRID/BIOGRID-ALL-3.4.158.tab2.txt';bp = process_file(fname)

