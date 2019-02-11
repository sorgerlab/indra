__all__ = ['process_pmc']

import json
import logging
import requests


logger = logging.getLogger(__name__)


default_output_fname = 'rlimsp_output.json'
RLIMSP_URL = 'https://research.bioinformatics.udel.edu/itextmine/api/data/rlims/pmc'


class RLIMSP_Error(Exception):
    pass


def process_pmc(pmc_id, output_fname=default_output_fname, with_grounding=True):
    """Get an output from RLIMS-p for the given pmic id."""
    if with_grounding:
        resp = requests.get(RLIMSP_URL + '.normed/pmcid/%s' % pmc_id)
    else:
        resp = requests.get(RLIMSP_URL + '/pmcid/%s' % pmc_id)

    if resp.status_code != 200:
        raise RLIMSP_Error("Bad status code: %d - %s"
                           % (resp.status_code, resp.reason))

    return
