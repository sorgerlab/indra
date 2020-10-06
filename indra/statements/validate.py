"""This module implements a number of functions that can be used to
validate INDRA Statements."""
import re
from indra.databases.identifiers import identifiers_mappings, \
    non_grounding, non_registry, identifiers_registry


class UnknownNamespace(ValueError):
    pass


class InvalidIdentifier(ValueError):
    pass


def validate_ns(db_ns):
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry or db_ns in non_registry \
            or db_ns in non_grounding:
        return True
    return False


def assert_valid_ns(db_ns):
    if not validate_ns(db_ns):
        raise UnknownNamespace(db_ns)


def validate_id(db_ns, db_id):
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry:
        if re.match(identifiers_registry[identifiers_ns]['pattern'], db_id):
            return True
        else:
            return False
    elif db_ns in non_registry or db_ns in non_grounding:
        return True
    else:
        return False


def assert_valid_id(db_ns, db_id):
    if not validate_id(db_ns, db_id):
        raise InvalidIdentifier(f'{db_ns}:{db_id}')


def assert_valid_db_refs(db_refs):
    for db_ns, db_id in db_refs.items():
        assert_valid_ns(db_ns)
        assert_valid_id(db_ns, db_id)


def validate_db_refs(db_refs):
    return all(validate_ns(db_ns) and validate_id(db_ns, db_id)
               for db_ns, db_id in db_refs.items())


def validate_statement(stmt):
    return all(validate_db_refs(agent.db_refs)
               for agent in stmt.real_agent_list())


def assert_valid_statement(stmt):
    for agent in stmt.real_agent_list():
        assert_valid_db_refs(agent.db_refs)
