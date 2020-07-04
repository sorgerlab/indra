import logging
from xml.etree import ElementTree
from .processor import DrugbankProcessor

logger = logging.getLogger(__name__)


def process_xml(fname):
    logger.info('Loading %s...' % fname)
    et = ElementTree.parse(fname)
    logger.info('Extracting DrugBank statements...')
    dp = DrugbankProcessor(et)
    dp.extract_inhibitions()
    return dp