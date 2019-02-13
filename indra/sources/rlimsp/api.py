__all__ = ['process_pmc']

import logging
import requests

from .processor import RlimspProcessor


logger = logging.getLogger(__name__)


RLIMSP_URL = 'https://research.bioinformatics.udel.edu/itextmine/api/data/rlims/pmc'


class RLIMSP_Error(Exception):
    pass


def process_pmc(pmcid, with_grounding=True):
    """Get an output from RLIMS-p for the given pmic id.

    Parameters
    ----------
    pmcid : str
        A PMCID, with the prefix PMC, of the paper to be "read".
    with_grounding : bool
        The RLIMS-P web service provides two endpoints, one pre-grounded, the
        other not so much. The grounded endpoint returns far less content, and
        may perform some grounding that can be handled by the grounding mapper.
    """
    if with_grounding:
        resp = requests.get(RLIMSP_URL + '.normed/pmcid/%s' % pmcid)
    else:
        resp = requests.get(RLIMSP_URL + '/pmcid/%s' % pmcid)

    if resp.status_code != 200:
        raise RLIMSP_Error("Bad status code: %d - %s"
                           % (resp.status_code, resp.reason))

    rp = RlimspProcessor(resp.json())

    return rp
