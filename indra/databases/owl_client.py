"""A client for OWL-sourced identifier mappings."""

import json
import os
import pickle
from collections import defaultdict
from operator import attrgetter
from typing import Any, Mapping, TYPE_CHECKING

from tqdm import tqdm

from indra.databases.obo_client import OntologyClient, prune_empty_entries
from indra.resources import get_resource_path

if TYPE_CHECKING:
    import pronto


class OwlClient(OntologyClient):
    """A base client for data that's been grabbed via OWL."""

    @staticmethod
    def entry_from_term(term: "pronto.Term", prefix,
                        remove_prefix: bool = False) -> Mapping[str, Any]:
        """Create a data dictionary from a Pronto term."""
        rels_dict = defaultdict(list)
        xrefs = []
        for xref in term.xrefs:
            try:
                xref_db, xref_id = xref.id.split(":", maxsplit=1)
            except ValueError:
                continue
            else:
                xrefs.append(dict(namespace=xref_db, id=xref_id))
        for child in term.subclasses(distance=1, with_self=False):
            child_db, child_id = child.id.split(':', maxsplit=1)
            if remove_prefix and child_db.lower() == prefix.lower():
                rels_dict["is_a"].append(child_id)
            else:
                rels_dict["is_a"].append(child.id)

        term_ns, term_id = term.id.split(':', maxsplit=1)
        term_ns = term_ns.lower()
        return {
            "namespace": term_ns,
            "id": term_id,
            "name": term.name,
            "synonyms": [s.description for s in term.synonyms],
            "xrefs": xrefs,
            "alt_ids": sorted(term.alternate_ids),
            "relations": dict(rels_dict),
        }

    @classmethod
    def entries_from_ontology(
        cls,
        prefix: str,
        ontology: "pronto.Ontology",
        *,
        skip_obsolete: bool = True,
        remove_prefix: bool = False,
    ):
        prefix = prefix.upper()
        rv = []
        for term in tqdm(ontology.terms(), desc=f"[{prefix}]"):
            if term.obsolete and skip_obsolete:
                continue
            if not term.id.startswith(prefix):
                continue
            rv.append(cls.entry_from_term(term, prefix,
                                          remove_prefix=remove_prefix))
        return rv

    @classmethod
    def update_resource(
        cls,
        prefix: str,
        ontology: "pronto.Ontology",
        skip_obsolete: bool = True,
        remove_prefix: bool = False,
    ):
        prefix = prefix.lower()
        entries = cls.entries_from_ontology(
            prefix=prefix, ontology=ontology, skip_obsolete=skip_obsolete,
            remove_prefix=remove_prefix
        )
        entries = prune_empty_entries(
            entries,
            {"synonyms", "xrefs", "alt_ids", "relations"},
        )
        entries = sorted(
            entries,
            key=attrgetter("id") if remove_prefix else _id_key,
        )

        resource_path = get_resource_path(f"{prefix}.json")
        with open(resource_path, "w") as file:
            json.dump(entries, file, indent=1, sort_keys=True)

    @classmethod
    def update_from_obo_library(
        cls,
        prefix: str,
        extension: str = "owl",
        **kwargs,
    ):
        prefix = prefix.lower()
        cache_path = get_resource_path(f"{prefix}.{extension}.pkl")

        if os.path.exists(cache_path):
            with open(cache_path, "rb") as file:
                ontology = pickle.load(file)
        else:
            try:
                import pronto
            except ImportError:
                raise ImportError(
                    "To use the INDRA OWL Client, you must first"
                    "install Pronto with `pip install pronto`."
                )
            ontology = pronto.Ontology.from_obo_library(
                f"{prefix.upper()}.{extension}")
            with open(cache_path, "wb") as file:
                pickle.dump(ontology, file, protocol=pickle.HIGHEST_PROTOCOL)

        cls.update_resource(prefix=prefix, ontology=ontology, **kwargs)

    @classmethod
    def update_from_file(
        cls,
        prefix: str,
        file,
        **kwargs,
    ):
        try:
            import pronto
        except ImportError:
            raise ImportError(
                "To use the INDRA OWL Client, you must first"
                "install Pronto with `pip install pronto`."
            )
        ontology = pronto.Ontology(file)
        cls.update_resource(prefix=prefix, ontology=ontology, **kwargs)


def _id_key(x):
    return int(x["id"].split(':')[1])


if __name__ == "__main__":
    OwlClient.update_from_obo_library("ido")
