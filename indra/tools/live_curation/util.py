import json
import logging
from indra.statements import Statement
from . import default_key_base

logger = logging.getLogger(__name__)


def _clean_key(s3key):
    # Check if default_key_base ('indra_models') is present in key
    s3key = s3key if default_key_base in s3key else \
        default_key_base + '/' + s3key

    # Replace double slashes
    s3key = s3key.replace('//', '/')

    # Ommit file part of key, assume it ends with json if it is present
    s3key = '/'.join([s for s in s3key.split('/')[:-1]]) if \
        s3key.endswith('.json') else s3key

    # Ensure last char in string is not '/'
    s3key = s3key[:-1] if s3key.endswith('/') else s3key

    return s3key


def _json_loader(fpath):
    logger.info('Loading json file %s' % fpath)
    with open(fpath, 'r') as f:
        return json.load(f)


def _json_dumper(jsonobj, fpath):
    try:
        logger.info('Saving json object to file %s' % fpath)
        with open(fpath, 'w') as f:
            json.dump(obj=jsonobj, fp=f, indent=1)
        return True
    except Exception as e:
        logger.error('Could not save json')
        logger.exception(e)
        return False


def _json_to_stmts_dict(stmt_jsons):
    """Return dict of statements keyed by uuid's from json statements

    This function is the inverse of _stmts_dict_to_json()

    Parameters
    ----------
    stmt_jsons : list(json)
        A list of json statements

    Returns
    -------
    dict
        Dict with statements keyed by their uuid's: {uuid: stmt}
    """
    loaded_stmts = [Statement._from_json(s) for s in stmt_jsons]
    return {s.uuid: s for s in loaded_stmts}


def _stmts_dict_to_json(stmt_dict):
    """Make a json representation from dict of statements

    This function is the inverse of _json_to_stmts_dict()

    Parameters
    ----------
    stmt_dict : dict
        Dict with statements keyed by their uuid's: {uuid: stmt}

    Returns
    -------
    list(json)
        A list of json statements
    """
    return [s.to_json() for _, s in stmt_dict.items()]


