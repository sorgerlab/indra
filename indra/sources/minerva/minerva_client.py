import json
import requests
import csv
from collections import defaultdict

default_map_name = 'covid19map'
base_url = 'https://%s.elixir-luxembourg.org/minerva/api/'
resource_url = ('https://git-r3lab.uni.lu/covid/models/-/raw/master/'
                'Integration/MINERVA_build/resources.csv')


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
        'bioEntities/elements/?columns=id,name,type,elementId,complexId,references'
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
    references = [(ref.get('type'), ref.get('resource')) for ref in refs]
    if element.get('name'):
        references.append(('TEXT', element['name']))
    return references


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


def get_ids_to_refs(model_id, map_name=default_map_name):
    config = get_config(map_name)
    project_id = get_project_id_from_config(config)
    model_elements = get_model_elements(model_id, project_id, map_name)
    object_ids_to_element_ids = {}
    ids_to_refs = {}
    complex_ids_to_members = defaultdict(set)
    for element in model_elements:
        object_ids_to_element_ids[element['id']] = element['elementId']
        ref = get_element_references(element)
        if ref:
            ids_to_refs[element['elementId']] = ref
        if element.get('complexId'):
            complex_ids_to_members[element['complexId']].add(
                element['elementId'])
    complex_members = {}
    for complex_id, members in complex_ids_to_members.items():
        complex_members[object_ids_to_element_ids[complex_id]] = members
    return ids_to_refs, complex_members


def get_model_ids(map_name=default_map_name):
    config = get_config(map_name)
    project_id = get_project_id_from_config(config)
    models = get_models(project_id, map_name)
    model_names_to_ids = {}
    for model in models:
        model_names_to_ids[model['name']] = model['idObject']
    return model_names_to_ids


def get_sif_filenames_to_ids(map_name=default_map_name):
    model_names_to_ids = get_model_ids()
    filenames_to_ids = {}
    res = requests.get(resource_url)
    csv_reader = csv.reader(res.text.splitlines())
    for row in csv_reader:
        model_name = row[3]
        if model_name in model_names_to_ids:
            fname = row[1].split('/')[-1][:-4] + '_raw.sif'
            filenames_to_ids[fname] = model_names_to_ids[model_name]
    return filenames_to_ids
