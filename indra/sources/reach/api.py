"""Methods for obtaining a reach processor containing indra statements.

Many file formats are supported. Many will run reach.
"""
import json
import logging
import requests

from indra.literature import id_lookup
import indra.literature.pmc_client as pmc_client
import indra.literature.pubmed_client as pubmed_client
from .processor import ReachProcessor


logger = logging.getLogger(__name__)

try:
    # For offline reading
    from .reader import ReachReader, ReachOfflineReadingError, JavaException
    reach_reader = ReachReader()
    try_offline = True
except Exception as e:
    logger.warning('Could not import jnius, offline reading option will not '
                   'be available.')
    logger.debug(e)
    try_offline = False

reach_text_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/text'
reach_nxml_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/nxml'
local_text_url = 'http://localhost:8080/api/text'
local_nxml_url = 'http://localhost:8080/api/uploadFile'
default_output_fname = 'reach_output.json'


def process_pmc(pmc_id, offline=False, url=None,
                output_fname=default_output_fname,
                organism_priority=None):
    """Return a ReachProcessor by processing a paper with a given PMC id.

    Uses the PMC client to obtain the full text. If it's not available,
    None is returned.

    Parameters
    ----------
    pmc_id : str
        The ID of a PubmedCentral article. The string may start with PMC but
        passing just the ID also works.
        Examples: 3717945, PMC3717945
        https://www.ncbi.nlm.nih.gov/pmc/
    offline : Optional[bool]
        If set to True, the REACH system is run offline via a JAR file.
        Otherwise (by default) the web service is called. Default: False
    url : Optional[str]
        URL for a REACH web service instance, which is used for reading if
        provided. If not provided but offline is set to False (its default
        value), the Arizona REACH web service is called
        (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
        Default: None
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    # Loading content from PMC first
    logger.info('Loading %s from PMC' % pmc_id)
    xml_str = pmc_client.get_xml(pmc_id)
    if xml_str is None:
        return None
    # Write into a file in the working folder
    fname = pmc_id + '.nxml'
    with open(fname, 'wb') as fh:
        fh.write(xml_str.encode('utf-8'))
    # Try to get the PMID for the paper so that the evidence pmid
    # attribute can be set correctly
    logger.info('Looking up PMID for %s' % pmc_id)
    ids = id_lookup(pmc_id, 'pmcid')
    pmid = ids.get('pmid')
    # Now process the NXML file with the provided arguments
    logger.info('Processing %s with REACH' % pmc_id)
    rp = process_nxml_file(fname, citation=pmid, offline=offline, url=url,
                           output_fname=output_fname,
                           organism_priority=organism_priority)
    return rp


def process_pubmed_abstract(pubmed_id, offline=False, url=None,
                            output_fname=default_output_fname, **kwargs):
    """Return a ReachProcessor by processing an abstract with a given Pubmed id.

    Uses the Pubmed client to get the abstract. If that fails, None is
    returned.

    Parameters
    ----------
    pubmed_id : str
        The ID of a Pubmed article. The string may start with PMID but
        passing just the ID also works.
        Examples: 27168024, PMID27168024
        https://www.ncbi.nlm.nih.gov/pubmed/
    offline : Optional[bool]
        If set to True, the REACH system is run offline via a JAR file.
        Otherwise (by default) the web service is called. Default: False
    url : Optional[str]
        URL for a REACH web service instance, which is used for reading if
        provided. If not provided but offline is set to False (its default
        value), the Arizona REACH web service is called
        (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
        Default: None
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.
    **kwargs : keyword arguments
        All other keyword arguments are passed directly to `process_text`.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    # Get the abstract from PubMed, if that fails, return None
    abs_txt = pubmed_client.get_abstract(pubmed_id)
    if abs_txt is None:
        return None
    # Process the text with the provided arguments
    rp = process_text(abs_txt, citation=pubmed_id, offline=offline, url=url,
                      output_fname=output_fname, **kwargs)
    # For some applications, the section type of the text is important so
    # that annotation is set here.
    if rp and rp.statements:
        for st in rp.statements:
            for ev in st.evidence:
                ev.epistemics['section_type'] = 'abstract'
    return rp


def process_text(text, citation=None, offline=False, url=None,
                 output_fname=default_output_fname, timeout=None,
                 organism_priority=None):
    """Return a ReachProcessor by processing the given text.

    Parameters
    ----------
    text : str
        The text to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. This is used when the text to be processed comes from
        a publication that is not otherwise identified. Default: None
    offline : Optional[bool]
        If set to True, the REACH system is run offline via a JAR file.
        Otherwise (by default) the web service is called. Default: False
    url : Optional[str]
        URL for a REACH web service instance, which is used for reading if
        provided. If not provided but offline is set to False (its default
        value), the Arizona REACH web service is called
        (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
        Default: None
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    timeout : Optional[float]
        This only applies when reading online (`offline=False`). Only wait for
        `timeout` seconds for the api to respond.
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        json_str = _read_content_offline(text, 'text')
    # If we are not reading offline then the old and new service interfaces
    # are the same so we can use a shared function
    else:
        if url is None:
            url = reach_text_url
        json_str = _read_text_service(text, url, timeout)

    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation=citation,
                                organism_priority=organism_priority)


def process_nxml_str(nxml_str, citation=None, offline=False,
                     url=None, output_fname=default_output_fname,
                     organism_priority=None):
    """Return a ReachProcessor by processing the given NXML string.

    NXML is the format used by PubmedCentral for papers in the open
    access subset.

    Parameters
    ----------
    nxml_str : str
        The NXML string to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. Default: None
    offline : Optional[bool]
        If set to True, the REACH system is run offline via a JAR file.
        Otherwise (by default) the web service is called. Default: False
    url : Optional[str]
        URL for a REACH web service instance, which is used for reading if
        provided. If not provided but offline is set to False (its default
        value), the Arizona REACH web service is called
        (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
        Default: None
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        json_str = _read_content_offline(nxml_str, 'nxml')
    else:
        # Use the Arizona URL by default if not given
        if url is None:
            url = reach_nxml_url
        # Print warning but proceed with reading
        if url == reach_nxml_url:
            logger.warning('Remote REACH webservice might get stuck when ' +
                           'processing NXML. Running local instance of REACH' +
                           ' is recommended.')
            json_str = _read_nxml_str_service_old(nxml_str, url)
        # Otherwise we assume that the web service is more recent than the
        # Arizona one and requires the new protocol.
        else:
            with open('temp_file.nxml', 'wb') as f:
                f.write(nxml_str.encode('utf-8'))
            json_str = _read_nxml_file_service_new('temp_file.nxml', url)

    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation=citation,
                                organism_priority=organism_priority)


def process_nxml_file(file_name, citation=None, offline=False,
                      url=None, output_fname=default_output_fname,
                      organism_priority=None):
    """Return a ReachProcessor by processing the given NXML file.

    NXML is the format used by PubmedCentral for papers in the open
    access subset.

    Parameters
    ----------
    file_name : str
        The name of the NXML file to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. Default: None
    offline : Optional[bool]
        If set to True, the REACH system is run offline via a JAR file.
        Otherwise (by default) the web service is called. Default: False
    url : Optional[str]
        URL for a REACH web service instance, which is used for reading if
        provided. If not provided but offline is set to False (its default
        value), the Arizona REACH web service is called
        (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
        Default: None
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    # First, if we are reading offline, we read the file and proceed
    if offline:
        with open(file_name, 'rb') as f:
            nxml_str = f.read().decode('utf-8')
            json_str = _read_content_offline(nxml_str, 'nxml')
    # If we are using the Arizona service, we use the old protocol
    elif url is None or url == reach_nxml_url:
        json_str = _read_nxml_file_service_old(file_name, url=reach_nxml_url)
    # Otherwise we use the new protocol
    else:
        json_str = _read_nxml_file_service_new(file_name, url=url)
    # Finally, we process the JSON output
    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation=citation,
                                organism_priority=organism_priority)


def process_json_file(file_name, citation=None, organism_priority=None):
    """Return a ReachProcessor by processing the given REACH json file.

    The output from the REACH parser is in this json format. This function is
    useful if the output is saved as a file and needs to be processed.
    For more information on the format, see: https://github.com/clulab/reach

    Parameters
    ----------
    file_name : str
        The name of the json file to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. Default: None
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_str(json_str, citation=citation,
                                    organism_priority=organism_priority)
    except IOError:
        logger.error('Could not read file %s.' % file_name)


def process_json_str(json_str, citation=None, organism_priority=None):
    """Return a ReachProcessor by processing the given REACH json string.

    The output from the REACH parser is in this json format.
    For more information on the format, see: https://github.com/clulab/reach

    Parameters
    ----------
    json_str : str
        The json string to be processed.
    citation : Optional[str]
        A PubMed ID passed to be used in the evidence for the extracted INDRA
        Statements. Default: None
    organism_priority : Optional[list of str]
        A list of Taxonomy IDs providing prioritization among organisms
        when choosing protein grounding. If not given, the default behavior
        takes the first match produced by Reach, which is prioritized to be
        a human protein if such a match exists.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    fields = ['frame-id', 'argument-label', 'object-meta',
              'doc-id', 'is-hypothesis', 'is-negated',
              'is-direct', 'found-by']
    for field in fields:
        json_str = json_str.replace(field, field.replace('-', '_'))
    try:
        json_dict = json.loads(json_str)
    except ValueError as e:
        logger.error('Could not decode JSON string.')
        logger.exception(e)
        return None
    rp = ReachProcessor(json_dict, pmid=citation,
                        organism_priority=organism_priority)
    rp.get_modifications()
    rp.get_complexes()
    rp.get_activation()
    rp.get_translocation()
    rp.get_regulate_amounts()
    rp.get_conversion()
    return rp


def _read_content_offline(content, content_type='text'):
    """Return a json string by processing the given text with offline
    REACH reader.

    Parameters
    ----------
    content : str
        The text to be processed.
    content_type : str
        Whether the content is a regular text or NXML.

    Returns
    -------
    json_str : bytes
        The json string produced by REACH reader.
    """
    if not try_offline:
        logger.error('Offline reading is not available.')
        return None
    try:
        api_ruler = reach_reader.get_api_ruler()
    except ReachOfflineReadingError as e:
        logger.error(e)
        logger.error('Cannot read offline because the REACH ApiRuler '
                     'could not be instantiated.')
        return None
    try:
        if content_type == 'text':
            result_map = api_ruler.annotateText(content, 'fries')
        elif content_type == 'nxml':
            result_map = api_ruler.annotateNxml(content, 'fries')
        else:
            raise ValueError('Invalid content_type: %s' % content_type)
    except JavaException as e:
        logger.error('Could not process %s.' % content_type)
        logger.error(e)
        return None
    # REACH version < 1.3.3
    json_str = result_map.get('resultJson')
    if not json_str:
        # REACH version >= 1.3.3
        json_str = result_map.get('result')
    if json_str is None:
        logger.warning('No results retrieved')
        return None
    if not isinstance(json_str, bytes):
        json_str = json_str.encode('utf-8')
    return json_str


def _read_text_service(text, url=reach_text_url, timeout=None):
    """Return a json string by processing the given text with online REACH API.

    Parameters
    ----------
    text : str
        The text to be processed.
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
    timeout : Optional[float]
        Only wait for `timeout` seconds for the api to respond.

    Returns
    -------
    json_str : bytes
        The json string returned by REACH API.
    """
    params = {'text': text.encode('utf-8')}
    try:
        res = requests.post(url, params=params, timeout=timeout)
    except requests.exceptions.RequestException as e:
        logger.error('Could not connect to REACH service:')
        logger.error(e)
        return None
    # TODO: we could use res.json() here to get a dict
    # directly
    # This is a byte string
    json_str = res.content
    return json_str


def _read_nxml_file_service_old(nxml_file, url=reach_nxml_url):
    with open(nxml_file, 'r', encoding='utf8') as fh:
        nxml_str = fh.read()
    return _read_nxml_str_service_old(nxml_str, url=url)


def _read_nxml_str_service_old(nxml_str, url=reach_nxml_url):
    """Return a json string by processing the given NXML string with remote
    REACH webservice.

    Parameters
    ----------
    nxml_str : str
        The NXML string to be processed.
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.

    Returns
    -------
    json_str : bytes
        The json string returned by REACH API.
    """
    data = {'nxml': nxml_str}
    try:
        res = requests.post(url, data)
    except requests.exceptions.RequestException as e:
        logger.error('Could not connect to REACH service:')
        logger.error(e)
        return None
    if res.status_code != 200:
        logger.error('Could not process NXML via REACH service.'
                     + 'Status code: %d' % res.status_code)
        return None
    json_str = res.content
    return json_str


def _read_nxml_file_service_new(file_name, url=local_nxml_url):
    """Return a json string by processing the given NXML file with locally
    running instance of REACH webservice.

    Parameters
    ----------
    file_name : str
        The name of the NXML file to be processed.
    url : Optional[str]
        URL for REACH service. By default, localhost on port 8080 is called.

    Returns
    -------
    json_str : bytes
        The json string returned by REACH API.
    """
    with open(file_name, 'rb') as f:
        try:
            res = requests.post(url, files={'file': f})
        except requests.exceptions.RequestException as e:
            logger.error('Could not connect to REACH service:')
            logger.error(e)
            return None
    if res.status_code != 200:
        logger.error('Could not process NXML via REACH service.'
                     + 'Status code: %d' % res.status_code)
        return None
    json_str = res.content
    return json_str
