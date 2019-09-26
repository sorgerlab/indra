from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from indra.sources.cwms.processor import CWMSProcessor
from indra.sources.cwms.rdf_processor import CWMSRDFProcessor
from indra.sources.trips import client

logger = logging.getLogger(__name__)


def process_text(text, save_xml='cwms_output.xml'):
    """Processes text using the CWMS web service.

    Parameters
    ----------
    text : str
        Text to process

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    xml = client.send_query(text, 'cwmsreader')

    # There are actually two EKBs in the xml document. Extract the second.
    first_end = xml.find('</ekb>')  # End of first EKB
    second_start = xml.find('<ekb', first_end)  # Start of second EKB
    second_end = xml.find('</ekb>', second_start)  # End of second EKB
    second_ekb = xml[second_start:second_end+len('</ekb>')]  # second EKB
    if save_xml:
        with open(save_xml, 'wb') as fh:
            fh.write(second_ekb.encode('utf-8'))
    return process_ekb(second_ekb)


def process_ekb_file(fname):
    """Processes an EKB file produced by CWMS.

    Parameters
    ----------
    fname : str
        Path to the EKB file to process.

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    # Process EKB XML file into statements
    with open(fname, 'rb') as fh:
        ekb_str = fh.read().decode('utf-8')
    return process_ekb(ekb_str)


def process_ekb(ekb_str):
    """Processes an EKB string produced by CWMS.

    Parameters
    ----------
    ekb_str : str
        EKB string to process

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    # Process EKB XML into statements
    cp = CWMSProcessor(ekb_str)
    cp.extract_causal_relations()
    cp.extract_correlations()
    cp.extract_migrations()
    cp.extract_events()
    return cp


def process_rdf_file(text, rdf_filename):
    """Process CWMS's RDF output for the given statement and returns a
    processor populated with INDRA statements.

    Parameters
    ----------
    text : str
        Sentence to process
    rdf_filename : str
        The RDF filename to process

    Returns
    -------
    cp : indra.sources.cwms.CWMSRDFProcessor
        A CWMSProcessor instance, which contains a list of INDRA Statements
        as its statements attribute.
    """
    cp = CWMSRDFProcessor(text, rdf_filename)
    return cp
