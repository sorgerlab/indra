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

        self.entries = {}
        self.alt_to_id = {}
        self.name_to_id = {}
        self.synonym_to_id = {}

        with open(self.mapping_path) as file:
            entries = json.load(file)

        self.entries = {entry['id']: entry for entry in entries}

        ambig_synonyms = set()
        for db_id, entry in self.entries.items():
            xrs = defaultdict(list)
            for xref in entry['xrefs']:
                xrs[xref['namespace']].append(xref['id'])
            entry['xrefs'] = dict(xrs)

            self.name_to_id[entry['name']] = db_id
            for synonym in entry['synonyms']:
                # Make a note of this is an ambiguous synonym so that we can
                # get rid of it after the loop, e.g., "multiciliation"
                if synonym in self.synonym_to_id:
                    ambig_synonyms.add(synonym)
                self.synonym_to_id[synonym] = db_id

            for db_alt_id in entry['alt_ids']:
                if db_alt_id in self.entries:
                    raise ValueError(
                        'Problem with integrity of {}:{}'.format(
                            self.prefix, db_alt_id
                        )
                    )
                self.alt_to_id[db_alt_id] = db_id
        # Remove all ambiguous synonyms
        self.synonym_to_id = {k: v for k, v in self.synonym_to_id.items()
                              if k not in ambig_synonyms}

    @staticmethod
    def entries_from_graph(obo_graph, prefix, remove_prefix=False,
                           allowed_synonyms=None):
        """Return processed entries from an OBO graph."""
        allowed_synonyms = allowed_synonyms if allowed_synonyms is not None \
            else {'EXACT', 'RELATED'}

        prefix_upper = prefix.upper()
        entries = []
        for node, data in obo_graph.nodes(data=True):
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

            # For simplicity, here we only take rels from the same ontology
            # but in principle, we could consider ones across ontologies
            rels_dict = defaultdict(list)
            if 'is_a' in data:
                rels_dict['is_a'] = data.get('is_a')
            for rel in data.get('relationship', []):
                rel_type, target = rel.split(' ', maxsplit=1)
                rels_dict[rel_type].append(target)
            for rel_type, rels in rels_dict.items():
                rel_own = [entry for entry in
                           sorted(set(rels)) if entry.startswith(prefix_upper)]
                rel_own = [(entry if not remove_prefix
                            else entry.split(':', maxsplit=1)[1])
                           for entry in rel_own]
                rels_dict[rel_type] = rel_own
            rels_dict = dict(rels_dict)

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
                'relations': rels_dict,
            })
        return entries

    @staticmethod
    def update_resource(directory, url, prefix, *args, remove_prefix=False,
                        allowed_synonyms=None):
        """Write the OBO information to files in the given directory."""
        resource_path = _make_resource_path(directory, prefix)
        obo_path = os.path.join(directory, '%s.obo.pkl' % prefix)
        if os.path.exists(obo_path):
            with open(obo_path, 'rb') as file:
                g = pickle.load(file)
        else:
            g = obonet.read_obo(url)
            with open(obo_path, 'wb') as file:
                pickle.dump(g, file)

        entries = \
            OboClient.entries_from_graph(g, prefix=prefix,
                                         remove_prefix=remove_prefix,
                                         allowed_synonyms=allowed_synonyms)
        entries = prune_empty_entries(entries,
                                      {'synonyms', 'xrefs',
                                       'alt_ids', 'relations'})
        with open(resource_path, 'w') as file:
            json.dump(entries, file, indent=1, sort_keys=True)

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
        db_name : str or None
            The name corresponding to the given ID.
        """
        return self.entries.get(db_id, {}).get('name')

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

    def get_id_from_name_or_synonym(self, txt):
        """Return the database id corresponding to the given name or synonym.

        Note that the way the OboClient is constructed, ambiguous synonyms are
        filtered out. Further, this function prioritizes names over synonyms
        (i.e., it first looks up the ID by name, and only if that fails,
        it attempts a synonym-based lookup). Overall, these mappings are
        guaranteed to be many-to-one.

        Parameters
        ----------
        txt : str
            The name or synonym to be converted.

        Returns
        -------
        db_id : str
            The ID corresponding to the given name or synonym.
        """
        name_id = self.get_id_from_name(txt)
        if name_id:
            return name_id
        return self.synonym_to_id.get(txt)

    def get_id_from_alt_id(self, db_alt_id):
        """Return the canonical database id corresponding to the alt id.

        Parameters
        ----------
        db_alt_id : str
            The alt id to be converted.

        Returns
        -------
        db_id : str or None
            The ID corresponding to the given alt id.
        """
        return self.alt_to_id.get(db_alt_id)

    def get_relations(self, db_id):
        """Return the isa relationships corresponding to a given ID.

        Parameters
        ----------
        db_id : str
            The ID whose isa relationships should be returned

        Returns
        -------
        dict
            A dict keyed by relation type with each entry a list of IDs of the
            terms that are in the given relation with the given ID.
        """
        return self.entries.get(db_id, {})

    def get_relation(self, db_id, rel_type):
        """Return the isa relationships corresponding to a given ID.

        Parameters
        ----------
        db_id : str
            The ID whose isa relationships should be returned
        rel_type : str
            The type of relationships to get, e.g., is_a, part_of

        Returns
        -------
        list of str
            The IDs of the terms that are in the given relation with the given
            ID.
        """
        return self.entries.get(db_id, {}).get(rel_type, [])


def prune_empty_entries(entries, keys):
    for entry in entries:
        for key in keys:
            if key in entry and not entry[key]:
                entry.pop(key)
    return entries