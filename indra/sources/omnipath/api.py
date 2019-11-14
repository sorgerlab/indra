import logging
import requests
from .processor import OmniPathLiganReceptorProcessor, \
    OmniPathModificationProcessor

logger = logging.getLogger("omnipath")


op_url = 'http://omnipathdb.org'


def process_from_web():
    ptm_json = _get_modifications()
    return OmniPathModificationProcessor(ptm_json)


def process_from_pypath(reload_resources=False, force=False):
    """Get all receptor ligand interactions from the omnipath pypath module

    Parameters
    ----------
    pa : pypath.main.PyPath
        An instance of a PyPath object containing the network
        representing ligand-receptor interactions

    Returns
    -------
    stmts : list[indra.statements.Statement]
        A list of indra statements"""
    # Import here rather than globally to remove overhead
    from pypath import main as pypath_main, data_formats
    pa = pypath_main.PyPath()
    pa.init_network(data_formats.ligand_receptor)

    if reload_resources:
        # Todo wipe the cache (stored in ~/.pypath/cache) clean and
        #  re-download the resources. Warn the user that it takes a lot of
        #  time to download it
        pass

    return OmniPathLiganReceptorProcessor(pa)


def _get_modifications():
    """Get all PTMs from Omnipath in JSON format.

    Returns
    -------
    JSON content for PTMs.
    """
    #params = {'format': 'json', 'fields':['sources', 'references']}
    params = {'format': 'json', 'fields':['sources']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    if not res.status_code == 200 or not res.text:
        return None
    else:
        return res.json()

