import os
import csv
from indra.statements.validate import validate_db_refs, validate_ns
from indra.preassembler.grounding_mapper import default_grounding_map
from indra.preassembler.grounding_mapper import default_misgrounding_map


def test_misgrounding_map_entries():
    bad_entries = []
    for text, db_refs in default_misgrounding_map.items():
        if not validate_db_refs(db_refs):
            bad_entries.append([text, db_refs])
    assert not bad_entries, bad_entries


def test_grounding_map_entries():
    bad_entries = []
    for text, db_refs in default_grounding_map.items():
        if not validate_db_refs(db_refs):
            bad_entries.append([text, db_refs])
    assert not bad_entries, bad_entries
