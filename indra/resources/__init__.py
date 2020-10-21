"""This module contains a number of resource files that INDRA uses
to perform tasks such as name standardization and ID mapping."""
import os
import json

RESOURCES_PATH = os.path.dirname(os.path.abspath(__file__))


def open_resource_file(resource_name, *args, **kwargs):
    """Return a file handle to an INDRA resource file."""
    if not resource_name.startswith(RESOURCES_PATH):
        resource_path = os.path.join(RESOURCES_PATH, resource_name)
    else:
        resource_path = resource_name

    return open(resource_path, *args, **kwargs)


def get_resource_path(fname):
    """Return the absolute path to a file in the resource folder."""
    return os.path.join(RESOURCES_PATH, fname)


def load_resource_json(fname):
    """Load a given JSON file from the resources folder.

    Parameters
    ----------
    fname: str
        The name of the json file in the resources folder.

    Returns
    -------
    json
        The content of the JSON file loaded into a dict/list.
    """
    with open(get_resource_path(fname), 'r') as fh:
        return json.load(fh)
