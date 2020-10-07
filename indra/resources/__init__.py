import os

RESOURCES_PATH = os.path.dirname(os.path.abspath(__file__))


def open_resource_file(resource_name, *args, **kwargs):
    """Return a file handle to an INDRA resource file."""
    if not resource_name.startswith(RESOURCES_PATH):
        resource_path = os.path.join(RESOURCES_PATH, resource_name)
    else:
        resource_path = resource_name

    return open(resource_path, *args, **kwargs)
