import requests
import logging
from .processor import SifProcessor
from .minerva_client import get_sif_filenames_to_ids


logger = logging.getLogger()
base_url = ('https://git-r3lab.uni.lu/covid/models/-/raw/master/'
            'Executable%20Modules/SBML_qual_build/sif/')


def process_file(filename, model_id):
    """Get statements by processing a single local SIF file.

    Parameters
    ----------
    filename : str
        A name (or path) of a local SIF file to process.
    model_id : int
        ID of a model corresponding to file content. Model ID is needed to
        find relevant references.

    Returns
    -------
    sp : indra.source.minerva.SifProcessor
        An instance of a SifProcessor with extracted INDRA statements.
    """
    return process_files({model_id: filename})


def process_files(ids_to_filenames):
    """Get statements by processing one or more local SIF files.

    Parameters
    ----------
    ids_to_file_names : dict
        A dictionary mapping model IDs to files containing model content as
        SIF. Model IDs are needed to find relevant references.

    Returns
    -------
    sp : indra.source.minerva.SifProcessor
        An instance of a SifProcessor with extracted INDRA statements.
    """
    model_id_to_sif_strs = {}
    for model_id, filename in ids_to_filenames.items():
        with open(filename, 'r') as f:
            sif_strs = f.readlines()
        model_id_to_sif_strs[model_id] = sif_strs
    return process_sif_strs(model_id_to_sif_strs)


def process_from_web(filenames='all'):
    """Get statements by processing remote SIF files.

    Parameters
    ----------
    filenames : list or str('all')
        Filenames for models that need to be processed (for full list of
        available models see
        https://git-r3lab.uni.lu/covid/models/-/tree/master/
        Executable%20Modules/SBML_qual_build/sif). If set to 'all'
        (default), then all available models will be processed.

    Returns
    -------
    sp : indra.source.minerva.SifProcessor
        An instance of a SifProcessor with extracted INDRA statements.
    """
    filenames_to_ids = get_sif_filenames_to_ids()
    if filenames == 'all':
        filenames = list(filenames_to_ids.keys())
    model_id_to_sif_strs = {}
    for fname in filenames:
        model_id = filenames_to_ids[fname]
        url = base_url + fname
        res = requests.get(url)
        if res.status_code == 200:
            sif_strs = res.text.split('\n')
            model_id_to_sif_strs[model_id] = sif_strs
        else:
            logger.warning('Could not get content from file %s, skipping '
                           'model %d' % (fname, model_id))
    return process_sif_strs(model_id_to_sif_strs)


def process_sif_strs(model_id_to_sif_strs):
    sp = SifProcessor(model_id_to_sif_strs)
    sp.extract_statements()
    return sp
