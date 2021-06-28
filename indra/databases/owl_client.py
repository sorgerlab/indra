"""A client for OWL-sourced identifier mappings."""

import json
import os
import pickle
from collections import defaultdict
from typing import Any, Mapping, TYPE_CHECKING

from tqdm import tqdm

from indra.databases.obo_client import OntologyClient, prune_empty_entries
from indra.resources import get_resource_path

if TYPE_CHECKING:
    import pronto


class OwlClient(OntologyClient):
    """A base client for data that's been grabbed via OWL."""

    @staticmethod
    def entry_from_term(term: "pronto.Term") -> Mapping[str, Any]:
        """Create a data dictionary from a Pronto term."""
        rels_dict = defaultdict(list)
        xrefs = []
        for xref in term.xrefs:
            try:
                xref_db, xref_id = xref.id.split(":")
            except ValueError:
                continue
            else:
                xrefs.append(dict(namespace=xref_db, id=xref_id))
        for child in term.subclasses(distance=1, with_self=False):
            rels_dict["is_a"].append(child.id)

        namespace, identifier = term.id.split(":")

        return {
            "namespace": namespace,
            "id": identifier,
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
    ):
        prefix = prefix.upper()
        rv = []
        for term in tqdm(ontology.terms(), desc=f"[{prefix}]"):
            if term.obsolete and skip_obsolete:
                continue
            if not term.id.startswith(prefix):
                continue
            rv.append(cls.entry_from_term(term))
        return rv

    @classmethod
    def update_resource(
        cls,
        prefix: str,
        ontology: "pronto.Ontology",
        skip_obsolete: bool = True,
    ):
        prefix = prefix.lower()
        entries = cls.entries_from_ontology(
            prefix=prefix, ontology=ontology, skip_obsolete=skip_obsolete
        )
        entries = prune_empty_entries(
            entries,
            {"synonyms", "xrefs", "alt_ids", "relations"},
        )
        entries = sorted(entries, key=lambda x: int(x["id"]))

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
            with cache_path.open("rb") as file:
                ontology = pickle.load(file)
        else:
            try:
                import pronto
            except ImportError:
                raise ImportError(
                    "To use the INDRA OWL Client, you must first"
                    "install Pronto with `pip install pronto`."
                )
            ontology = pronto.Ontology.from_obo_library(f"{prefix.upper()}.{extension}")
            with cache_path.open("wb") as file:
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


if __name__ == "__main__":
    OwlClient.update_from_obo_library("ido")
