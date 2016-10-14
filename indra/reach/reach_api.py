from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import json
import logging
import tempfile
import requests
import indra.literature.pmc_client as pmc_client
import indra.literature.pubmed_client as pubmed_client
from indra.reach.processor import ReachProcessor
# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('reach')

try:
    # For offline reading
    from indra.java_vm import autoclass, JavaException
    from indra.reach.reach_reader import ReachReader
    reach_reader = ReachReader()
    try_offline = True
except Exception:
    logger.error('Could not import jnius, offline reading cannot be used.')
    try_offline = False

reach_text_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/text'
reach_nxml_url = 'http://agathon.sista.arizona.edu:8080/odinweb/api/nxml'


def process_pmc(pmc_id, offline=False):
    """Return a ReachProcessor by processing a paper with a given PMC id.

    Uses the PMC client to obtain the full text. If it's not available,
    None is returned.

    Parameters
    ----------
    pmc_id : str
        The ID of a PubmedCentral article. The string may start with PMC but
        passing just the ID also works.
        Examples: 3717945, PMC3717945
        http://www.ncbi.nlm.nih.gov/pmc/
    offline : Optional[bool]
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    xml_str = pmc_client.get_xml(pmc_id)
    if xml_str is None:
        return None
    fname = pmc_id + '.nxml'
    with open(fname, 'wb') as fh:
        fh.write(xml_str.encode('utf-8'))
    rp = process_nxml_file(fname, citation=pmc_id, offline=offline)
    return rp


def process_pubmed_abstract(pubmed_id, offline=False):
    """Return a ReachProcessor by processing an abstract with a given Pubmed id.

    Uses the Pubmed client to get the abstract. If that fails, None is
    returned.

    Parameters
    ----------
    pubmed_id : str
        The ID of a Pubmed article. The string may start with PMID but
        passing just the ID also works.
        Examples: 27168024, PMID27168024
        http://www.ncbi.nlm.nih.gov/pubmed/
    offline : Optional[bool]
        If set to True, the REACH system is ran offline. Otherwise (by default)
        the web service is called. Default: False

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    abs_txt = pubmed_client.get_abstract(pubmed_id)
    if abs_txt is None:
        return None
    rp = process_text(abs_txt, citation=pubmed_id, offline=offline)
    return rp


def process_text(text, citation=None, offline=False):
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

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        if not try_offline:
            logger.error('Offline reading is not available.')
            return None
        api_ruler = reach_reader.get_api_ruler()
        if api_ruler is None:
            logger.error('Cannot read offline because the REACH ApiRuler ' + \
                         'could not be instantiated.')
            return None
        try:
            result_map = api_ruler.annotateText(text, 'fries')
        except JavaException as e:
            logger.error('Could not process text.')
            logger.error(e)
            return None
        json_str = result_map.get('resultJson')
    else:
        data = {'text': text.encode('utf-8')}
        try:
            res = requests.post(reach_text_url, data)
        except requests.exceptions.RequestException as e:
            logger.error('Could not connect to REACH service:')
            logger.error(e)
            return None
        # TODO: we could use res.json() here to get a dict 
        # directly
        # This is a unicode string
        json_str = res.text
        # In Python2 this is a test against future str type
        assert isinstance(json_str, str)
    with open('reach_output.json', 'wt') as fh:
        fh.write(json_str)
    return process_json_str(json_str, citation)


def process_nxml_str(nxml_str, citation=None, offline=False):
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

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    if offline:
        fname = 'tmp.nxml'
        with open(fname, 'wb') as fh:
            fh.write(nxml_str.encode('utf-8'))
        rp = process_nxml_file(fname, citation, True)
        return rp
    else:
        data = {'nxml': nxml_str}
        try:
            res = requests.post(reach_nxml_url, data)
        except requests.exceptions.RequestException as e:
            logger.error('Could not connect to REACH service:')
            logger.error(e)
            return None
        if res.status_code != 200:
            logger.error('Could not process NXML via REACH service.' + \
                         'Status code: %d' % res.status_code)
            return None
        json_str = res.text
        with open('reach_output.json', 'wt') as fh:
            fh.write(json_str)
        return process_json_str(json_str, citation)


def process_nxml_file(file_name, citation=None, offline=False):
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

    Returns
    -------
    rp : ReachProcessor
        A ReachProcessor containing the extracted INDRA Statements
        in rp.statements.
    """
    # Offline we use the API ruler directly to read the nxml fle
    if offline:
        if not try_offline:
            logger.error('Offline reading is not available.')
            return None
        api_ruler = reach_reader.get_api_ruler()
        if api_ruler is None:
            logger.error('Cannot read offline because the REACH ApiRuler' +\
                         'could not be instantiated.')
            return None
        try:
            result_map = api_ruler.annotateNxml(file_name, 'fries')
        except JavaException as e:
            logger.error('Could not process NXML.')
            logger.error(e)
            return None
        json_str = result_map.get('resultJson')
        return process_json_str(json_str)
    # For the web service, we read the file and process it as a string
    else:
        with open(file_name, 'rb') as f:
            nxml_str = f.read().decode('utf-8')
        return process_nxml_str(nxml_str, citation, False)


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
        with open(file_name, 'rt') as fh:
            json_str = fh.read()
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
    assert isinstance(json_str, basestring)
    json_str = json_str.replace('frame-id','frame_id')
    json_str = json_str.replace('argument-label','argument_label')
    json_str = json_str.replace('object-meta','object_meta')
    json_str = json_str.replace('doc-id','doc_id')
    json_str = json_str.replace('is-hypothesis','is_hypothesis')
    json_str = json_str.replace('is-negated','is_negated')
    json_str = json_str.replace('is-direct','is_direct')
    json_str = json_str.replace('found-by','found_by')
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
    return rp
