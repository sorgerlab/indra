__all__ = ['process_from_webservice', 'process_from_json_file',
           'process_from_jsonish_str']

import json
import logging
import requests

from .processor import RlimspProcessor


logger = logging.getLogger(__name__)


RLIMSP_URL = ('https://research.bioinformatics.udel.edu/itextmine/api/data/'
              'rlims/')


class RLIMSP_Error(Exception):
    pass


def process_from_webservice(id_val, id_type='pmcid', source='pmc'):
    """Return an output from RLIMS-p for the given PubMed ID or PMC ID.

    The web service is documented at: https://research.bioinformatics.udel.edu/itextmine/api/.
    The /data/rlims URL endpoint is extended with three additional elements:
    /{collection}/{key}/{value} where collection is "medline" or "pmc", key is
    "pmid" or "pmcid", and value is a specific PMID or PMCID.

    Parameters
    ----------
    id_val : str
        A PMCID, with the prefix PMC, or PMID, with no prefix, of the paper to
        be "read". Corresponds to the "value" argument of the REST API.
    id_type : Optional[str]
        Either 'pmid' or 'pmcid'. The default is 'pmcid'. Corresponds to the
        "key" argument of the REST API.
    source : Optional[str]
        Either 'pmc' or 'medline', whether you want pmc fulltext or medline
        abstracts. Corresponds to the "collection" argument of the REST API.

    Returns
    -------
    :py:class:`indra.sources.rlimsp.processor.RlimspProcessor`
        An RlimspProcessor which contains a list of extracted INDRA Statements
        in its statements attribute.
    """
    resp = requests.get(RLIMSP_URL + '%s/%s/%s' % (source, id_type, id_val))

    if resp.status_code != 200:
        raise RLIMSP_Error("Bad status code: %d - %s"
                           % (resp.status_code, resp.reason))

    rp = RlimspProcessor(resp.json())
    rp.extract_statements()
    return rp


def process_from_json_file(filename, doc_id_type=None):
    """Process RLIMSP extractions from a bulk-download JSON file.

    Parameters
    ----------
    filename : str
        Path to the JSON file.
    doc_id_type : Optional[str]
        In some cases the RLIMS-P paragraph info doesn't contain 'pmid' or
        'pmcid' explicitly, instead if contains a 'docId' key. This parameter
        allows defining what ID type 'docId' sould be interpreted as. Its
        values should be 'pmid' or 'pmcid' or None if not used.

    Returns
    -------
    :py:class:`indra.sources.rlimsp.processor.RlimspProcessor`
        An RlimspProcessor which contains a list of extracted INDRA Statements
        in its statements attribute.
    """
    with open(filename, 'rt') as f:
        lines = f.readlines()
        json_list = []
        for line in lines:
            json_list.append(json.loads(line))
        rp = RlimspProcessor(json_list, doc_id_type=doc_id_type)
        rp.extract_statements()
    return rp


def process_from_jsonish_str(jsonish_str, doc_id_type=None):
    """Process RLIMSP extractions from a bulk-download JSON file.

    Parameters
    ----------
    jsonish_str : str
        The contents of one of the not-quite-json files you can find here:
        https://hershey.dbi.udel.edu/textmining/export
    doc_id_type : Optional[str]
        In some cases the RLIMS-P paragraph info doesn't contain 'pmid' or
        'pmcid' explicitly, instead if contains a 'docId' key. This parameter
        allows defining what ID type 'docId' sould be interpreted as. Its
        values should be 'pmid' or 'pmcid' or None if not used.

    Returns
    -------
    :py:class:`indra.sources.rlimsp.processor.RlimspProcessor`
        An RlimspProcessor which contains a list of extracted INDRA Statements
        in its statements attribute.
    """
    lines = jsonish_str.splitlines()
    json_list = []
    for line in lines:
        json_list.append(json.loads(line))
    rp = RlimspProcessor(json_list, doc_id_type=doc_id_type)
    rp.extract_statements()
    return rp
