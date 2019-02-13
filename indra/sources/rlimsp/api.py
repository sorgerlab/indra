__all__ = ['process_pmc']

import logging
import requests

from .processor import RlimspProcessor


logger = logging.getLogger(__name__)


RLIMSP_URL = 'https://research.bioinformatics.udel.edu/itextmine/api/data/rlims/pmc'


class RLIMSP_Error(Exception):
    pass


def process_pmc(pmc_id, with_grounding=True):
    """Get an output from RLIMS-p for the given pmic id."""
    if with_grounding:
        resp = requests.get(RLIMSP_URL + '.normed/pmcid/%s' % pmc_id)
    else:
        resp = requests.get(RLIMSP_URL + '/pmcid/%s' % pmc_id)

    if resp.status_code != 200:
        raise RLIMSP_Error("Bad status code: %d - %s"
                           % (resp.status_code, resp.reason))

    rp = RlimspProcessor(resp.json())

    return rp
