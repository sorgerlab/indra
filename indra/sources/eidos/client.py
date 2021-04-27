__all__ = ['process_text']

import requests


def process_text(text, webservice):
    """Process a given text with an Eidos webservice at the given address.

    Note that in most cases this function should not be used directly, rather,
    used indirectly by calling `indra.sources.eidos.process_text` with
    the webservice parameter.

    Parameters
    ----------
    text : str
        The text to be read using Eidos.
    webservice : str
        The address where the Eidos web service is running, e.g.,
        http://localhost:9000.

    Returns
    -------
    dict
        A JSON dict of the results from the Eidos webservice.
    """
    res = requests.post('%s/process_text' % webservice,
                        json={'text': text})
    res.raise_for_status()
    json_dict = res.json()
    return json_dict
