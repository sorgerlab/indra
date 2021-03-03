import json
import requests

default_map_name = 'covid19map'
base_url = 'https://%s.elixir-luxembourg.org/minerva/api/'


def get_config(map_name=default_map_name):
    url = (base_url % map_name) + 'configuration/'
    res = requests.get(url)
    res.raise_for_status()
    return res.json()


def get_project_id_from_config(config):
    options = config.get('options', [])
    for option in options:
        if option.get('type') == 'DEFAULT_MAP':
            return option.get('value')
    return None


def get_models(project_id, map_name=default_map_name):
    url = (base_url % map_name) + ('projects/%s/models/' % project_id)
    res = requests.get(url)
    res.raise_for_status()
    return res.json()


def get_model_elements(model_id, project_id, map_name=default_map_name):
    url = (base_url % map_name) + \
        ('projects/%s/models/%s/' % (project_id, model_id)) + \
        'bioEntities/elements/?columns=id,name,type,references'
    res = requests.get(url)
    res.raise_for_status()
    return res.json()


def get_all_model_elements(models, project_id, map_name=default_map_name):
    all_elements = []
    for model in models:
        model_id = model['idObject']
        model_elements = get_model_elements(model_id, project_id,
                                            map_name)
        all_elements += model_elements
    return all_elements


def get_element_references(element):
    refs = element.get('references', [])
    return [(ref.get('type'), ref.get('resource')) for ref in refs]


def get_all_valid_element_refs(map_name=default_map_name):
    config = get_config(map_name)
    project_id = get_project_id_from_config(config)
    models = get_models(project_id, map_name)
    all_model_elements = get_all_model_elements(models, project_id,
                                                map_name)
    element_refs = [get_element_references(element) for element
                    in all_model_elements]
    valid_element_refs = [ref for ref in element_refs if ref]
    return valid_element_refs


