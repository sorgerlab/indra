import logging
from xml.etree import ElementTree
from .processor import DrugbankProcessor

logger = logging.getLogger(__name__)


def process_xml(fname):
    """Return a processor by extracting Statement from DrugBank XML.

    Parameters
    ----------
    fname : str
        The path to a DrugBank XML file to process.

    Returns
    -------
    DrugbankProcessor
        A DrugbankProcessor instance which contains a list of INDRA
        Statements in its statements attribute that were extracted
        from the given XML file.
    """
    logger.info('Loading %s...' % fname)
    et = ElementTree.parse(fname)
    logger.info('Extracting DrugBank statements...')
    dp = DrugbankProcessor(et)
    dp.extract_inhibitions()
    return dp