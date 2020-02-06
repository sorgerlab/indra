__all__ = ['process_text', 'reground']

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


def reground(ont_yml, texts, webservice, topk=10, is_canonicalized=False,
             filter=True):
    """Ground concept texts given an ontology with an Eidos web service.

    Parameters
    ----------
    ont_yml : str
        A serialized YAML string representing the ontology.
    texts : list[str]
        A list of concept texts to ground.
    webservice : str
        The address where the Eidos web service is running, e.g.,
        http://localhost:9000.
    topk : Optional[int]
        The number of top scoring groundings to return. Default: 10
    is_canonicalized : Optional[bool]
        If True, the texts are assumed to be canonicalized. If False,
        Eidos will canonicalize the texts which yields much better groundings
        but is slower. Default: False
    filter : Optional[bool]
        If True, Eidos filters the ontology to remove determiners from examples
        and other similar operations. Should typically be set to True.
        Default: True

    Returns
    -------
    dict
        A JSON dict of the results from the Eidos webservice.
    """
    params = {
        'ontologyYaml': ont_yml,
        'texts': texts,
        'topk': topk,
        'isAlreadyCanonicalized': is_canonicalized,
        'filter': filter
    }
    res = requests.post('%s/reground' % webservice,
                        json=params)
    res.raise_for_status()
    json_dict = res.json()
    return json_dict
