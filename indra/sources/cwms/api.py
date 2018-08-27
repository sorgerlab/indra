import logging
from indra.sources.cwms.processor import CWMSProcessor
from indra.sources.cwms.rdf_processor import CWMSRDFProcessor
from indra.sources.trips import trips_client

logger = logging.getLogger('cwms')


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
    xml = trips_client.send_query(text, 'cwmsreader')

    # There are actually two EKBs in the xml document. Extract the second.
    first_end = xml.find('</ekb>')  # End of first EKB
    second_start = xml.find('<ekb', first_end)  # Start of second EKB
    second_end = xml.find('</ekb>', second_start)  # End of second EKB
    second_ekb = xml[second_start:second_end+len('</ekb>')]  # second EKB
    if save_xml:
        with open(save_xml, 'wb') as fh:
            fh.write(second_ekb.encode('utf-8'))
    return process_ekb(second_ekb)


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
