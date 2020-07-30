import json
import logging

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
