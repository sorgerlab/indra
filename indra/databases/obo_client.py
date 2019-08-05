"""A client for OBO-sourced identifier mappings."""

import os
import pickle

import obonet

__all__ = [
    'OboClient',
    'RESOURCES',
]

HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(HERE, os.pardir, 'resources')


def _make_resource_path(directory, prefix):
    return os.path.join(directory, '{prefix}.tsv'.format(prefix=prefix))


class OboClient:
    """A base client for data that's been grabbed via OBO"""

    def __init__(self, prefix, *, directory: str = RESOURCES):
        """Read the OBO file export at the given path."""
        self.prefix = prefix
        self.directory = directory
        self.mapping_path = _make_resource_path(self.directory, self.prefix)
        self.id_to_name = {}
        self.name_to_id = {}
        with open(self.mapping_path) as file:
            next(file)  # throw away the header
            for line in file:
                db_id, db_name = line.strip().split('\t')
                self.id_to_name[db_id] = db_name
                self.name_to_id[db_name] = db_id

    @staticmethod
    def update_resource(
            directory,
            url,
            prefix,
            *,
            remove_prefix: bool = False,
    ) -> None:
        """Write the OBO information to files in the given directory."""
        prefix_upper = prefix.upper()

        tsv_path = _make_resource_path(directory, prefix)
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

        with open(tsv_path, 'w') as file:
            print(
                '{prefix}_id'.format(prefix=prefix),
                'name',
                sep='\t',
                file=file,
            )
            for node, data in g.nodes(data=True):
                if node.startswith(prefix_upper):
                    if remove_prefix:
                        node = node[len(prefix) + 1:]
                    print(node, data['name'], sep='\t', file=file)

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
