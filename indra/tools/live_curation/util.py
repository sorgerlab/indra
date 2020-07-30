import json
import logging
from indra.statements import stmts_from_json, stmts_to_json

logger = logging.getLogger(__name__)


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
        Dict with statements keyed by their uuids: {uuid: stmt}
    """
    loaded_stmts = stmts_from_json(stmt_jsons)
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
    return stmts_to_json(list(stmt_dict.values()))
