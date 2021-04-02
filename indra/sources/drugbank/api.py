import logging
from typing import Optional, Sequence, Union
from xml.etree import ElementTree
from .processor import DrugbankProcessor

logger = logging.getLogger(__name__)


def process_from_web(
    username: Optional[str] = None,
    password: Optional[str] = None,
    version: Optional[str] = None,
    prefix: Union[None, str, Sequence[str]] = None,
) -> DrugbankProcessor:
    """Get a processor using :func:`process_xml` with :mod:`drugbank_downloader`.

    Parameters
    ----------
    username :
        The DrugBank username. If not passed, looks up in the environment
        ``DRUGBANK_USERNAME``. If not found, raises a ValueError.
    password :
        The DrugBank password. If not passed, looks up in the environment
        ``DRUGBANK_PASSWORD``. If not found, raises a ValueError.
    version :
        The DrugBank version. If not passed, uses :mod:`bioversions` to
        look up the most recent version.
    prefix :
        The prefix and subkeys passed to :func:`pystow.ensure` to specify
        a non-default location to download the data to.

    Returns
    -------
    DrugbankProcessor
        A DrugbankProcessor instance which contains a list of INDRA
        Statements in its statements attribute that were extracted
        from the given DrugBank version
    """
    import drugbank_downloader
    et = drugbank_downloader.parse_drugbank(
        username=username,
        password=password,
        version=version,
        prefix=prefix,
    )
    return process_element_tree(et)


def process_xml(fname):
    """Return a processor by extracting Statements from DrugBank XML.

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
    return process_element_tree(et)


def process_element_tree(et):
    """Return a processor by extracting Statement from DrugBank XML.

    Parameters
    ----------
    et : xml.etree.ElementTree
        An ElementTree loaded from the DrugBank XML file to process.

    Returns
    -------
    DrugbankProcessor
        A DrugbankProcessor instance which contains a list of INDRA
        Statements in its statements attribute that were extracted
        from the given ElementTree.
    """
    logger.info('Extracting DrugBank statements...')
    dp = DrugbankProcessor(et)
    dp.extract_statements()
    return dp
