"""A client for OBO-sourced identifier mappings."""

import json
import logging
import os
import pickle
import re
from collections import Counter, defaultdict

import obonet

__all__ = [
    'OboClient',
    'RESOURCES',
]

HERE = os.path.dirname(os.path.abspath(__file__))
RESOURCES = os.path.join(HERE, os.pardir, 'resources')

logger = logging.getLogger(__name__)


def _make_resource_path(directory, prefix):
    return os.path.join(directory, '{prefix}.json'.format(prefix=prefix))


class OboClient:
    """A base client for data that's been grabbed via OBO"""

    def __init__(self, prefix, *, directory=RESOURCES):
        """Read the OBO file export at the given path."""
        self.prefix = prefix
        self.directory = directory
        self.mapping_path = _make_resource_path(self.directory, self.prefix)

        self.id_to_name = {}
        self.alt_to_id = {}
        self.name_to_id = {}
        self.id_to_xrefs = defaultdict(lambda: defaultdict(list))
        self.id_to_isa = {}

        with open(self.mapping_path) as file:
            entries = json.load(file)

        for entry in entries:
            db_id, db_name = entry['id'], entry['name']
            self.id_to_name[db_id] = db_name
            self.name_to_id[db_name] = db_id
            for xref in entry['xrefs']:
                xref_db, xref_db_id = xref['namespace'], xref['id']
                self.id_to_xrefs[db_id][xref_db].append(xref_db_id)

            for db_alt_id in entry['alt_ids']:
                if db_alt_id in self.id_to_name:
                    raise ValueError(
                        'Problem with integrity of {}:{}'.format(
                            self.prefix, db_alt_id
                        )
                    )
                self.alt_to_id[db_alt_id] = db_id
            self.id_to_isa[db_id] = entry['is_a']

    @staticmethod
    def update_resource(directory, url, prefix, *args, remove_prefix=False,
                        allowed_synonyms=None):
        """Write the OBO information to files in the given directory."""
        allowed_synonyms = allowed_synonyms if allowed_synonyms is not None \
            else {'EXACT', 'RELATED'}
        prefix_upper = prefix.upper()

        resource_path = _make_resource_path(directory, prefix)
        obo_path = os.path.join(directory, '%s.obo.pkl' % prefix)
        if os.path.exists(obo_path):
            with open(obo_path, 'rb') as file:
                g = pickle.load(file)
        else:
            g = obonet.read_obo(url)
            with open(obo_path, 'wb') as file:
                pickle.dump(g, file)

        entries = []
        for node, data in g.nodes(data=True):
            # There are entries in some OBOs that are actually from other
            # ontologies
            if not node.startswith(prefix_upper):
                continue
            if remove_prefix:
                node = node[len(prefix) + 1:]

            xrefs = []
            for xref in data.get('xref', []):
                db, db_id = xref.split(':', maxsplit=1)
                # Example: for EFO, we have xrefs like
                # PERSON: James Malone
                db_id = db_id.lstrip()
                # Example: for HP, we have xrefs like
                # MEDDRA:10050185 "Palmoplantar pustulosis"
                if ' ' in db_id:
                    db_id = db_id.split()[0]
                    logging.debug(
                        'Likely labeled %s:%s xref: %s. Recovered %s:%s',
                        prefix, node, xref, db, db_id,
                    )

                xrefs.append(dict(namespace=db, id=db_id))

            # For simplicity, here we only take isa from the same ontology
            # but in principle, we could consider ones across ontologies
            isa_raw = data.get('is_a', [])
            isa_own = [entry for entry in
                       sorted(set(isa_raw)) if entry.startswith(prefix_upper)]
            isa_own = [(entry if not remove_prefix
                        else entry.split(':', maxsplit=1)[1])
                       for entry in isa_own]

            synonyms = []
            for synonym in data.get('synonym', []):
                match = re.match(r'^\"(.+)\" (EXACT|RELATED|NARROW|BROAD)',
                                 synonym)
                syn, status = match.groups()
                if status in allowed_synonyms:
                    synonyms.append(syn)

            namespace = data.get('namespace', prefix)

            entries.append({
                'namespace': namespace,
                'id': node,
                'name': data['name'],
                'synonyms': synonyms,
                'xrefs': xrefs,
                'alt_ids': data.get('alt_id', []),
                'is_a': isa_own,
            })

        with open(resource_path, 'w') as file:
            json.dump(entries, file, indent=2, sort_keys=True)

    def count_xrefs(self):
        """Count how many xrefs there are to each database."""
        return Counter(
            xref_db
            for db_id, xref_map in self.id_to_xrefs.items()
            for xref_db, xref_db_ids in xref_map.items()
            for _ in xref_db_ids
        )

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

    def get_id_from_alt_id(self, db_alt_id):
        """Return the canonical database id corresponding to the alt id.

        Parameters
        ----------
        db_alt_id : str
            The alt id to be converted.

        Returns
        -------
        db_id : str
            The ID corresponding to the given alt id.
        """
        return self.alt_to_id.get(db_alt_id)

    def get_isa(self, db_id):
        return self.id_to_isa.get(db_id, [])
