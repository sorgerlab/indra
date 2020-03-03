"""Methods for obtaining a reach processor containing indra statements.

Many file formats are supported. Many will run reach.
"""
import json
import logging
import os
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


def process_pmc(pmc_id, offline=False, url=reach_nxml_url,
                output_fname=default_output_fname):
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
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
        This only applies when reading online (`offline=False`).
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    logger.info('Loading %s from PMC' % pmc_id)
    xml_str = pmc_client.get_xml(pmc_id)
    if xml_str is None:
        return None
    fname = pmc_id + '.nxml'
    with open(fname, 'wb') as fh:
        fh.write(xml_str.encode('utf-8'))
    logger.info('Looking up PMID for %s' % pmc_id)
    ids = id_lookup(pmc_id, 'pmcid')
    pmid = ids.get('pmid')
    logger.info('Processing %s with REACH' % pmc_id)
    rp = process_nxml_file(fname, citation=pmid, offline=offline, url=url,
                           output_fname=output_fname)
    return rp


def process_pubmed_abstract(pubmed_id, offline=False, url=reach_nxml_url,
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
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
        This only applies when reading online (`offline=False`).
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    **kwargs : keyword arguments
        All other keyword arguments are passed directly to `process_text`.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    abs_txt = pubmed_client.get_abstract(pubmed_id)
    if abs_txt is None:
        return None
    rp = process_text(abs_txt, citation=pubmed_id, offline=offline, url=url,
                      output_fname=output_fname, **kwargs)
    if rp and rp.statements:
        for st in rp.statements:
            for ev in st.evidence:
                ev.epistemics['section_type'] = 'abstract'
    return rp


def process_text(text, citation=None, offline=False, url=reach_text_url,
                 output_fname=default_output_fname, timeout=None):
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
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
        This only applies when reading online (`offline=False`).
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.
    timeout : Optional[float]
        This only applies when reading online (`offline=False`). Only wait for
        `timeout` seconds for the api to respond.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        json_str = get_json_str_offline(text, 'text')
    else:
        json_str = text_to_json_str_online(text, url, timeout)

    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation)


def process_nxml_str(nxml_str, citation=None, offline=False,
                     url=reach_nxml_url, output_fname=default_output_fname):
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
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
        This only applies when reading online (`offline=False`).
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        json_str = get_json_str_offline(nxml_str, 'nxml')
    else:
        if url == reach_nxml_url:
            json_str = nxml_to_json_str_remote(nxml_str, url)
        else:
            with open('temp_file.nxml', 'wb') as f:
                f.write(nxml_str)
            json_str = nxml_file_to_json_str_local('temp_file.nxml', url)

    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation)


def process_nxml_file(file_name, citation=None, offline=False,
                      url=reach_nxml_url, output_fname=default_output_fname):
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
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False
    url : Optional[str]
        URL for REACH service. By default, Arizona REACH web service is called.
        This only applies when reading online (`offline=False`).
    output_fname : Optional[str]
        The file to output the REACH JSON output to.
        Defaults to reach_output.json in current working directory.

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline or url == reach_nxml_url:
        with open(file_name, 'rb') as f:
            nxml_str = f.read().decode('utf-8')
            return process_nxml_str(nxml_str, citation, offline=offline,
                                    url=url, output_fname=output_fname)

    json_str = nxml_file_to_json_str_local(file_name, url)
    if json_str:
        with open(output_fname, 'wb') as fh:
            fh.write(json_str)
        return process_json_str(json_str.decode('utf-8'), citation)


def process_json_file(file_name, citation=None):
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

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    try:
        with open(file_name, 'rb') as fh:
            json_str = fh.read().decode('utf-8')
            return process_json_str(json_str, citation)
    except IOError:
        logger.error('Could not read file %s.' % file_name)


def process_json_str(json_str, citation=None):
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

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if not isinstance(json_str, str):
        raise TypeError('{} is {} instead of {}'.format(json_str,
                                                        json_str.__class__,
                                                        str))

    json_str = json_str.replace('frame-id', 'frame_id')
    json_str = json_str.replace('argument-label', 'argument_label')
    json_str = json_str.replace('object-meta', 'object_meta')
    json_str = json_str.replace('doc-id', 'doc_id')
    json_str = json_str.replace('is-hypothesis', 'is_hypothesis')
    json_str = json_str.replace('is-negated', 'is_negated')
    json_str = json_str.replace('is-direct', 'is_direct')
    json_str = json_str.replace('found-by', 'found_by')
    try:
        json_dict = json.loads(json_str)
    except ValueError:
        logger.error('Could not decode JSON string.')
        return None
    rp = ReachProcessor(json_dict, citation)
    rp.get_modifications()
    rp.get_complexes()
    rp.get_activation()
    rp.get_translocation()
    rp.get_regulate_amounts()
    return rp


def get_json_str_offline(content, content_type='text'):
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
    except JavaException as e:
        logger.error('Could not process %d.' % content_type)
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


def text_to_json_str_online(text, url=reach_text_url, timeout=None):
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
    if not isinstance(json_str, bytes):
        raise TypeError('{} is {} instead of {}'.format(
            json_str, json_str.__class__, bytes))
    return json_str


def nxml_to_json_str_remote(nxml_str, url=reach_nxml_url):
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
    json_str = res.text
    return res.text


def nxml_file_to_json_str_local(file_name, url=local_nxml_url):
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
    json_str = res.text
    return res.text
