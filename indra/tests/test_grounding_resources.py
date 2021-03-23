import os
import csv
from indra.statements.validate import validate_db_refs, validate_ns
from indra.preassembler.grounding_mapper import default_grounding_map
from indra.preassembler.grounding_mapper import default_misgrounding_map


# Namespaces that are not currently handled but still appear in statements
exceptions = ['CLO']

def test_misgrounding_map_entries():
    bad_entries = []
    for text, db_refs in default_misgrounding_map.items():
        if not validate_db_refs(db_refs):
            bad_entries.append([text, db_refs])
    assert not bad_entries, bad_entries


def test_grounding_map_entries():
    bad_entries = []
    for text, db_refs in default_grounding_map.items():
        if (not validate_db_refs(db_refs) and
            not (set(exceptions) & db_refs.keys())):
            bad_entries.append([text, db_refs])
    assert not bad_entries, bad_entries


def test_exceptional_unhandled():
    """Test that exceptional namespaces actually aren't handled.

    This will catch if we make an update that makes an exceptional namespace
    become a handled namespace. That way we can update the tests.
    """
    actually_handled = []
    for ns in exceptions:
        if validate_ns(ns):
            actually_handled.append(ns)
    assert not actually_handled, actually_handled






