import os
import logging
import requests
from .processor import OmniPathLiganReceptorProcessor, \
    OmniPathModificationProcessor

logger = logging.getLogger("omnipath")

try:
    from pypath import main as pypath_main, data_formats
    has_pypath = True
except ImportError:
    logger.info('PyPath is not available')
    pypath_main, data_formats = None, None
    has_pypath = False


op_url = 'http://omnipathdb.org'


def process_from_web():
    ptm_json = _get_modifications()
    return OmniPathModificationProcessor(ptm_json)


def process_from_pypath(reload_resources=False, force=False):
    """Get all receptor ligand interactions from the omnipath pypath module

    Parameters
    ----------
    reload_resources : bool
        If True, wipe the local cache (typically in ~/.pypath/cache),
        triggering a re-download of the resources.
    force : bool
        If True, don't ask user for permission to wipe the cache.

    Returns
    -------
    stmts : list[indra.statements.Statement]
        A list of indra statements"""
    if not has_pypath:
        logger.warning('Unable to run OmniPathLiganReceptorProcessor when '
                       'PyPath is not available')
        return None
    pa = pypath_main.PyPath()
    pa.init_network(data_formats.ligand_receptor)

    if reload_resources:
        success = _delete_omnipath_cache(force)
        if success:
            logger.info('Successfully emptied omnipath cache')
        else:
            logger.warning('Failed to empty cache')
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


def _delete_omnipath_cache(force=False):
    if not has_pypath:
        logger.warning('PyPath cache is not available: PyPath could not be '
                       'imported')
        return False
    from pypath.cache import get_cachedir
    cache_path = get_cachedir()
    if os.path.isdir(cache_path) and \
            len(os.walk(cache_path).__next__()[2]) > 0:
        logger.warning('Deleting the omnipath cache')
        if not force:
            print('Re-loading the omnipath resources can take up to an hour '
                  'for some of its resources.')
        ok = input('This action will remove all files in the omnipath '
                   'cahce. Proceed? [Y/n] ') if not force else 'n'
        try:
            if force or ok.lower() == 'y':
                for file in os.walk(cache_path).__next__()[2]:
                    file_path = os.path.join(cache_path, file)
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
        except Exception as e:
            logger.exception('Failed to delete file(s)')
            # Should raise the exception here, because if we partially
            # emptied the cache, we don't know if the needed resources are
            # there or not
            raise e
    else:
        logger.info('No files detected in %s' % cache_path)
        return False
