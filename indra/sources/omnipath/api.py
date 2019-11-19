import os
import logging
import requests
from .processor import OmniPathProcessor

logger = logging.getLogger("omnipath")


try:
    from pypath import main as pypath_main, data_formats
    from pypath.intera import Complex as pp_Complex
    has_pypath = True
except ImportError:
    logger.info('PyPath is not available')
    pypath_main = None
    data_formats = None
    pp_Complex = None
    has_pypath = False


op_url = 'http://omnipathdb.org'


def process_from_web():
    ptm_json = _get_modifications()
    return OmniPathProcessor(ptm_json=ptm_json)


def process_from_pypath(resources='all', reload_resources=False,
                        force=False):
    """Get all receptor ligand interactions from the omnipath pypath module

    Parameters
    ----------
    reload_resources : bool
        If True, wipe the local cache (typically in ~/.pypath/cache),
        triggering a re-download of the required resources.
    force : bool
        If True, don't ask user for permission to wipe the cache.

    Returns
    -------
    op : OmniPathProcessor
        An instance of OmniPathProcessor with statements in op.statements
    """
    if not has_pypath:
        logger.warning('Unable to process from PyPath: PyPath is not '
                       'available')
        return None

    if reload_resources:
        success = _delete_omnipath_cache(force)
        if success:
            logger.info('Successfully emptied omnipath cache')
        else:
            logger.warning('Failed to empty cache')

    ligrec_json = _get_ligand_receptor_interactions(resources)
    return OmniPathProcessor(ligrec_json=ligrec_json)


def _get_ligand_receptor_interactions(resources='all'):
    """

    Parameters
    ----------
    resources : set|'all'|str
        Either a string or a set of strings naming the resources to use.
        If any string is not part of the allowed resources, it will be
        ignored.

    Returns
    -------

    """
    allowed_resources = {'ramilowski2015', 'kirouac2010', 'hpmr',
                         'cellphonedb', 'guide2pharma'}

    # Helpers
    def _get_pypath_text_refs(article_id_list):
        text_refs = {}
        for ref in article_id_list:
            name = ref['idtype'].upper()
            try:
                id = int(ref['value'])
            except ValueError:
                id = ref['value']
            text_refs[name] = id
        return text_refs

    def _get_annotations(ref_info):
        """ref_info is a dict returned by the method 'info' of a
        pypath.refs.Reference object"""
        annotations = {}
        uid = ref_info['uids'][0]
        ref_dict = ref_info[uid]
        if ref_dict.get('recordstatus'):
            annotations['recordstatus'] = ref_dict['recordstatus']
        if ref_dict.get('pubstatus'):
            annotations['pubstatus'] = ref_dict['pubstatus']
        if ref_dict.get('pmcrefcount'):
            annotations['pmcrefcount'] = ref_dict['pmcrefcount']
        if ref_dict.get('reportnumber'):
            annotations['reportnumber'] = ref_dict['reportnumber']
        if ref_dict.get('availablefromurl'):
            annotations['availablefromurl'] = ref_dict['availablefromurl']
        if ref_dict.get('locationlabel'):
            annotations['locationlabel'] = ref_dict['locationlabel']
        return annotations

    def _get_complex_agents(up_id_string):
        if 'COMPLEX' in up_id_string:
            if ' ' in up_id_string:
                up_id_string = up_id_string.split()[-1]
            return [up_id for up_id in
                    up_id_string.split('COMPLEX:')[1].split('_')]
        return [up_id_string]

    def _get_participating_agents(pp_obj, s, t):
        """Return a list of UP ids of all the participating agents"""

        # Get participating agents
        if isinstance(pp_obj.vs[s]['name'], pp_Complex):
            # Catch the odd pypath.intera.Complex objects
            src_string = str(pp_obj.vs[s]['name'])
        else:
            src_string = pp_obj.vs[s]['name']
        source_agents = _get_complex_agents(src_string)

        if isinstance(pp_obj.vs[t]['name'], pp_Complex):
            # Catch the odd pypath.intera.Complex objects
            trg_string = str(pp_obj.vs[t]['name'])
        else:
            trg_string = pp_obj.vs[t]['name']
        target_agents = _get_complex_agents(trg_string)

        # Make unique list
        return list({*source_agents, *target_agents})

    # Main code body
    pa = pypath_main.PyPath()
    if isinstance(resources, str):
        if resources.lower() == 'all':
            pa.init_network(data_formats.ligand_receptor)
        elif resources in allowed_resources:
            pa.init_network({
                resources: data_formats.ligand_receptor[resources]
            })
        else:
            logger.warning('Resource %s not among allowed resources')
            return []
    elif isinstance(resources, (set, list)):
            resource_dict = {r: data_formats.ligand_receptor[r] for r in
                             resources if r in allowed_resources}
            pa.init_network(resource_dict)
    else:
        logger.error('resources not recognized as string or set of strings')
        return []

    ligrec_json = []

    for s, t in pa.graph.get_edgelist():
        agent_list = _get_participating_agents(pa, s, t)
        edge_obj = pa.get_edge(s, t)

        # Loop support
        for ref_name, ref_set in edge_obj['refs_by_source'].items():
            for ref_obj in ref_set:
                # Ref obj is a pypath.refs.Reference object
                # Check for PMID
                if ref_obj.pmid:
                    pmid = ref_obj.pmid
                else:
                    pmid = None
                try:
                    # pypath.refs.Reference.info() returns a dictionary
                    # of reference info. If load fails silently, an
                    # empty dict is returned.
                    ref_info = ref_obj.info()
                    if ref_info.get('uids'):
                        uid = ref_info['uids'][0]
                        # Text refs
                        text_refs = _get_pypath_text_refs(
                            ref_info[uid]['articleids'])
                        text_refs['nlmuniqueid'] = \
                            ref_info[uid]['nlmuniqueid']
                        text_refs['ISSN'] = ref_info[uid]['issn']
                        text_refs['ESSN'] = ref_info[uid]['essn']
                    else:
                        text_refs = None
                except TypeError as e:
                    logger.warning('Failed to load info')
                    logger.exception(e)
                    ref_info = None
                    text_refs = None

                # If both pmid and text_refs is None, skip this Complex
                if pmid is None and text_refs is None:
                    logger.info('Skipping statement for %s' % ref_name)
                    continue

                # Get annotations
                if ref_info:
                    annotations = _get_annotations(ref_info)
                    annotations['source_sub_id'] = ref_name
                else:
                    annotations = {'source_sub_id': ref_name}
                if ref_name == 'Ramilowski2015' and \
                        edge_obj['ramilowski_sources']:
                    annotations['ramilowski_sources'] = \
                        edge_obj['ramilowski_sources']
                if ref_name.lower() == 'cellphonedb' and \
                        edge_obj['cellphonedb_type']:
                    annotations['cellphonedb_type'] = \
                        edge_obj['cellphonedb_type']
                lr_entry = {
                    'pmid': pmid,
                    'annotations': annotations,
                    'text_refs': text_refs,
                    'agent_list': agent_list
                }
                ligrec_json.append(lr_entry)
    return ligrec_json


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
        logger.warning('PyPath cache is not available: PyPath could not'
                       ' be imported')
        return False
    from pypath.cache import get_cachedir
    cache_path = get_cachedir()
    if os.path.isdir(cache_path) and \
            len(os.walk(cache_path).__next__()[2]) > 0:
        logger.warning('Deleting the omnipath cache')
        if not force:
            print('Re-loading the omnipath resources can take up to an'
                  ' hour for some of its resources.')
        ok = input('This action will remove all files in the omnipath '
                   'cache. Proceed? [Y/n] ') if not force else 'n'
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
