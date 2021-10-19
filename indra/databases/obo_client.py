"""A client for OBO-sourced identifier mappings."""

import json
import logging
import os
import pathlib
import pickle
import re
from collections import Counter, defaultdict
from typing import List, Mapping, Optional

import obonet

from indra.resources import get_resource_path, load_resource_json

__all__ = [
    'OntologyClient',
    'OboClient',
    "PyOboClient",
]

HERE = pathlib.Path(__file__).parent.resolve()

logger = logging.getLogger(__name__)


class OntologyClient:
    """A base client class for OBO and OWL ontologies."""

    def __init__(self, prefix: str):
        """Read the OBO file export at the given path."""
        self.prefix = prefix.lower()
        self.entries = {
            entry['id']: entry for entry
            in load_resource_json(f'{prefix}.json')
        }
        self.alt_to_id = {}
        self.name_to_id = {}
        self.synonym_to_id = {}

        ambig_synonyms = set()
        for db_id, entry in self.entries.items():
            xrs = defaultdict(list)
            for xref in entry.get('xrefs', []):
                xrs[xref['namespace']].append(xref['id'])
            entry['xrefs'] = dict(xrs)

            self.name_to_id[entry['name']] = db_id
            for synonym in entry.get('synonyms', []):
                # Make a note of this is an ambiguous synonym so that we can
                # get rid of it after the loop, e.g., "multiciliation"
                if synonym in self.synonym_to_id:
                    ambig_synonyms.add(synonym)
                self.synonym_to_id[synonym] = db_id

            for db_alt_id in entry.get('alt_ids', []):
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

    def get_name_from_id(self, db_id: str) -> Optional[str]:
        """Return the database name corresponding to the given database ID.

        Parameters
        ----------
        db_id :
            The ID to be converted.

        Returns
        -------
        :
            The name corresponding to the given ID.
        """
        return self.entries.get(db_id, {}).get('name')

    def get_id_from_name(self, db_name: str) -> Optional[str]:
        """Return the database identifier corresponding to the given name.

        Parameters
        ----------
        db_name :
            The name to be converted.

        Returns
        -------
        :
            The ID corresponding to the given name.
        """
        return self.name_to_id.get(db_name)

    def get_id_from_name_or_synonym(self, txt: str) -> Optional[str]:
        """Return the database id corresponding to the given name or synonym.

        Note that the way the OboClient is constructed, ambiguous synonyms are
        filtered out. Further, this function prioritizes names over synonyms
        (i.e., it first looks up the ID by name, and only if that fails,
        it attempts a synonym-based lookup). Overall, these mappings are
        guaranteed to be many-to-one.

        Parameters
        ----------
        txt :
            The name or synonym to be converted.

        Returns
        -------
        :
            The ID corresponding to the given name or synonym.
        """
        name_id = self.get_id_from_name(txt)
        if name_id:
            return name_id
        return self.synonym_to_id.get(txt)

    def get_id_from_alt_id(self, db_alt_id: str) -> Optional[str]:
        """Return the canonical database id corresponding to the alt id.

        Parameters
        ----------
        db_alt_id :
            The alt id to be converted.

        Returns
        -------
        :
            The ID corresponding to the given alt id.
        """
        return self.alt_to_id.get(db_alt_id)

    def get_relations(self, db_id: str) -> Mapping[str, List[str]]:
        """Return the isa relationships corresponding to a given ID.

        Parameters
        ----------
        db_id :
            The ID whose isa relationships should be returned

        Returns
        -------
        :
            A dict keyed by relation type with each entry a list of IDs of the
            terms that are in the given relation with the given ID.
        """
        return self.entries.get(db_id, {})

    def get_relation(self, db_id: str, rel_type: str) -> List[str]:
        """Return the isa relationships corresponding to a given ID.

        Parameters
        ----------
        db_id :
            The ID whose isa relationships should be returned
        rel_type :
            The type of relationships to get, e.g., is_a, part_of

        Returns
        -------
        :
            The IDs of the terms that are in the given relation with the given
            ID.
        """
        return self.entries.get(db_id, {}).get(rel_type, [])


class OboClient(OntologyClient):
    """A base client for data that's been grabbed via OBO"""

    @staticmethod
    def entries_from_graph(obo_graph, prefix, remove_prefix=False,
                           allowed_synonyms=None, allowed_external_ns=None):
        """Return processed entries from an OBO graph."""
        allowed_synonyms = allowed_synonyms if allowed_synonyms is not None \
            else {'EXACT', 'RELATED'}

        prefix_upper = prefix.upper()
        entries = []

        for node, data in obo_graph.nodes(data=True):
            if 'name' not in data:
                continue
            # There are entries in some OBOs that are actually from other
            # ontologies. We either skip these entirely or if allowed
            # external name spaces are provided, we allow nodes that are
            # in one of those namespaces
            external_node = False
            if not node.startswith(prefix_upper):
                if allowed_external_ns and \
                        node.split(':')[0] in allowed_external_ns:
                    external_node = True
                else:
                    continue

            if not external_node and remove_prefix:
                node = node[len(prefix) + 1:]

            xrefs = []
            for xref in data.get('xref', []):
                try:
                    db, db_id = xref.split(':', maxsplit=1)
                # This is typically the case when the xref doesn't have
                # a separate name space in which case we skip it
                except ValueError:
                    continue
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
                           sorted(set(rels)) if entry.startswith(prefix_upper)
                           or (allowed_external_ns and
                               entry.split(':')[0] in allowed_external_ns)]
                rel_own = [(entry if ((not remove_prefix)
                                      or (allowed_external_ns
                                          and entry.split(':')[0] in
                                          allowed_external_ns))
                            else entry.split(':', maxsplit=1)[1])
                           for entry in rel_own]
                rels_dict[rel_type] = rel_own
            rels_dict = dict(rels_dict)

            synonyms = []
            for synonym in data.get('synonym', []):
                match = re.match(r'^\"(.+)\" (EXACT|RELATED|NARROW|BROAD|\[\])',
                                 synonym)
                syn, status = match.groups()
                if status == '[]':
                    status = 'EXACT'
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

    @classmethod
    def update_resource(cls, directory, url, prefix, *args, remove_prefix=False,
                        allowed_synonyms=None, allowed_external_ns=None):
        """Write the OBO information to files in the given directory."""
        resource_path = get_resource_path(f'{prefix}.json')
        obo_path = os.path.join(directory, '%s.obo.pkl' % prefix)
        if os.path.exists(obo_path):
            with open(obo_path, 'rb') as file:
                g = pickle.load(file)
        else:
            g = obonet.read_obo(url)
            with open(obo_path, 'wb') as file:
                pickle.dump(g, file)

        entries = \
            OboClient.entries_from_graph(
                g, prefix=prefix,
                remove_prefix=remove_prefix,
                allowed_synonyms=allowed_synonyms,
                allowed_external_ns=allowed_external_ns)
        entries = prune_empty_entries(entries,
                                      {'synonyms', 'xrefs',
                                       'alt_ids', 'relations'})

        def sort_key(x):
            val = x['id']
            if not remove_prefix:
                val = val.split(':')[1]
            try:
                val = int(val)
            except ValueError:
                pass
            return val

        entries = sorted(entries, key=sort_key)
        with open(resource_path, 'w') as file:
            json.dump(entries, file, indent=1, sort_keys=True)


class PyOboClient(OntologyClient):
    """A base client for data that's been grabbed via PyOBO."""

    @classmethod
    def update_by_prefix(cls, prefix: str):
        """Update the JSON data by looking up the ontology through PyOBO."""
        raise NotImplementedError


def prune_empty_entries(entries, keys):
    for entry in entries:
        for key in keys:
            if key in entry and not entry[key]:
                entry.pop(key)
    return entries
