__all__ = ['process_text', 'reground']

import requests


def process_text(text, webservice):
    res = requests.post('%s/process_text' % webservice,
                        json={'text': text})
    res.raise_for_status()
    json_dict = res.json()
    return json_dict


def reground(ont_yml, texts, webservice, topk=10, is_canonicalized=False,
             filter=True):
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
