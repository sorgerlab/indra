"""A client for OBO-sourced identifier mappings."""

import json
import os
import pickle
import re

import obonet

__all__ = [
    'OboClient',
    'RESOURCES',
]

HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(HERE, os.pardir, 'resources')


def _make_resource_path(directory, prefix):
    return os.path.join(directory, '{prefix}.json'.format(prefix=prefix))


OBO_SYNONYM = re.compile('EXACT|RELATED')


class OboClient:
    """A base client for data that's been grabbed via OBO"""

    def __init__(self, prefix, *, directory=RESOURCES):
        """Read the OBO file export at the given path."""
        self.prefix = prefix
        self.directory = directory
        self.mapping_path = _make_resource_path(self.directory, self.prefix)
        self.id_to_name = {}
        self.name_to_id = {}
        with open(self.mapping_path) as file:
            entries = json.load(file)

        for entry in entries:
            db_id, db_name = entry['id'], entry['name']
            self.id_to_name[db_id] = db_name
            self.name_to_id[db_name] = db_id

    @staticmethod
    def update_resource(directory, url, prefix, *args, remove_prefix=False):
        """Write the OBO information to files in the given directory."""
        prefix_upper = prefix.upper()

        resource_path = _make_resource_path(directory, prefix)
        obo_path = os.path.join(
            directory,
            '{prefix}.obo.pickle'.format(prefix=prefix),
        )
        if os.path.exists(obo_path):
            with open(obo_path, 'rb') as file:
                g = pickle.load(file)
        else:
            g = obonet.read_obo(url)
            with open(obo_path, 'wb') as file:
                pickle.dump(g, file)

        entries = []
        for node, data in g.nodes(data=True):
            if not node.startswith(prefix_upper):
                continue
            if remove_prefix:
                node = node[len(prefix) + 1:]
            entries.append({
                'namespace': prefix,
                'id': node,
                'name': data['name'],
                'synonyms': [
                    OBO_SYNONYM.split(synonym)[0].strip('" ')
                    for synonym in data.get('synonym', [])
                ],
                'xrefs': data.get('xref', []),
            })

        with open(resource_path, 'w') as file:
            json.dump(entries, file, indent=2, sort_keys=True)

    def get_name_from_id(self, db_id):
        """Return the database name corresponding to the given database ID.

        Parameters
        ----------
        db_id : str
            The ID to be converted.

        Returns
        -------
        db_name : str
            The name corresponding to the given ID.
        """
        return self.id_to_name.get(db_id)

    def get_id_from_name(self, db_name):
        """Return the database identifier corresponding to the given name.

        Parameters
        ----------
        db_name : str
            The name to be converted.

        Returns
        -------
        db_id : str
            The ID corresponding to the given name.
        """
        return self.name_to_id.get(db_name)
