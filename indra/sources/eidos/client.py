__all__ = ['process_text', 'reground_texts']

import tqdm
import requests
from indra.util import batch_iter


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


def reground_texts(texts, ont_yml, webservice, topk=10, is_canonicalized=False,
                   filter=True):
    """Ground concept texts given an ontology with an Eidos web service.

    Parameters
    ----------
    texts : list[str]
        A list of concept texts to ground.
    ont_yml : str
        A serialized YAML string representing the ontology.
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
    all_results = []
    for text_batch in tqdm.tqdm(batch_iter(texts, batch_size=500,
                                           return_func=list)):
        params = {
            'ontologyYaml': ont_yml,
            'texts': text_batch,
            'topk': topk,
            'isAlreadyCanonicalized': is_canonicalized,
            'filter': filter
        }
        res = requests.post('%s/reground' % webservice,
                            json=params)
        res.raise_for_status()
        all_results += grounding_dict_to_list(res.json())
    return all_results


def grounding_dict_to_list(groundings):
    """Transform the webservice response into a flat list."""
    all_grounding_lists = []
    for entry in groundings:
        grounding_list = []
        for grounding_dict in entry:
            grounding_list.append((grounding_dict['grounding'],
                                   grounding_dict['score']))
        grounding_list = sorted(grounding_list, key=lambda x: x[1],
                                reverse=True)
        all_grounding_lists.append(grounding_list)
    return all_grounding_lists
