import logging
from indra.sources.cwms.processor import CWMSProcessor, \
    CWMSProcessorCompositional
from indra.sources.trips import client

logger = logging.getLogger(__name__)

default_grounding_mode = 'flat'


def process_text(text, save_xml='cwms_output.xml',
                 extract_filter=None, grounding_mode=default_grounding_mode):
    """Processes text using the CWMS web service.

    Parameters
    ----------
    text : str
        Text to process
    save_xml : Optional[str]
        A file name in which to dump the output from CWMS.
        Default: cwms_output.xml
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence', 'association', 'event' and 'migration'.
        If not given, only Influences are extracted since processing other
        relation types can be time consuming. This argument can be used if
        the extraction of other relation types such as Events are also of
        interest.
    grounding_mode : Optional[str]
        Selects whether 'flat' or 'compositional' groundings should be
        extracted. Default: 'flat'.

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
    return process_ekb(second_ekb, extract_filter=extract_filter,
                       grounding_mode=grounding_mode)


def process_ekb_file(fname, extract_filter=None,
                     grounding_mode=default_grounding_mode):
    """Processes an EKB file produced by CWMS.

    Parameters
    ----------
    fname : str
        Path to the EKB file to process.
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence', 'association', 'event' and 'migration'.
        If not given, only Influences are extracted since processing other
        relation types can be time consuming. This argument can be used if
        the extraction of other relation types such as Events are also of
        interest.
    grounding_mode : Optional[str]
        Selects whether 'flat' or 'compositional' groundings should be
        extracted. Default: 'flat'.

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    # Process EKB XML file into statements
    with open(fname, 'rb') as fh:
        ekb_str = fh.read().decode('utf-8')
    return process_ekb(ekb_str, extract_filter=extract_filter,
                       grounding_mode=grounding_mode)


def process_ekb(ekb_str, extract_filter=None,
                grounding_mode=default_grounding_mode):
    """Processes an EKB string produced by CWMS.

    Parameters
    ----------
    ekb_str : str
        EKB string to process
    extract_filter : Optional[list]
        A list of relation types to extract. Valid values in the list are
        'influence', 'association', 'event' and 'migration'.
        If not given, only Influences are extracted since processing other
        relation types can be time consuming. This argument can be used if
        the extraction of other relation types such as Events are also of
        interest.
    grounding_mode : Optional[str]
        Selects whether 'flat' or 'compositional' groundings should be
        extracted. Default: 'flat'.

    Returns
    -------
    cp : indra.sources.cwms.CWMSProcessor
        A CWMSProcessor, which contains a list of INDRA statements in its
        statements attribute.
    """
    # Process EKB XML into statements
    if grounding_mode == 'flat':
        cp = CWMSProcessor(ekb_str)
    elif grounding_mode == 'compositional':
        cp = CWMSProcessorCompositional(ekb_str)
    else:
        raise ValueError('Invalid grounding mode: %s' % grounding_mode)
    if extract_filter is None or 'influence' in extract_filter:
        cp.extract_causal_relations()
    if extract_filter is not None and 'association' in extract_filter:
        cp.extract_correlations()
    if extract_filter is not None and 'migration' in extract_filter:
        cp.extract_migrations()
    if extract_filter is not None and 'event' in extract_filter:
        cp.extract_events()
    return cp
