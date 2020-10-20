"""This module implements a number of functions that can be used to
validate INDRA Statements."""
import re
from indra.statements import BioContext, RefContext, WorldContext
from indra.databases.identifiers import identifiers_mappings, \
    non_grounding, non_registry, identifiers_registry

text_ref_patterns = {
    'PMID': re.compile(r'(\d+)'),
    'PMCID': re.compile(r'PMC(\d+)'),
    # See https://www.crossref.org/blog/dois-and-matching-regular-expressions/
    # here I added a-z to allow non-capitalized DOI text, we could
    # change this is we want strict capital letters in the DOI
    'DOI':  re.compile(r'^10.\d{4,9}/[-._;()/:A-Za-z0-9]+$'),
}


class UnknownNamespace(ValueError):
    pass


class InvalidIdentifier(ValueError):
    pass


class InvalidTextRefs(ValueError):
    pass


class InvalidContext(ValueError):
    pass


def validate_ns(db_ns):
    """Return True if the given namespace is known.

    Parameters
    ----------
    db_ns : str
        The namespace.

    Returns
    -------
    bool
        True if the given namepsace is known, otherwise False.
    """
    try:
        assert_valid_ns(db_ns)
        return True
    except ValueError:
        return False


def assert_valid_ns(db_ns):
    """Raise UnknownNamespace error if the given namespace is unknown.

    Parameters
    ----------
    db_ns : str
        The namespace.
    """
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry or db_ns in non_registry \
            or db_ns in non_grounding:
        return
    raise UnknownNamespace(db_ns)


def validate_id(db_ns, db_id):
    """Return True if the given ID is valid in the given namespace.

    Parameters
    ----------
    db_ns : str
        The namespace.
    db_id : str
        The ID.

    Returns
    -------
    bool
        True if the given ID is valid in the given namespace.
    """
    try:
        assert_valid_id(db_ns, db_id)
        return True
    except ValueError:
        return False


def assert_valid_id(db_ns, db_id):
    """Raise InvalidIdentifier error if the ID is invalid in the given
    namespace.

    Parameters
    ----------
    db_ns : str
        The namespace.
    db_id : str
        The ID.
    """
    identifiers_ns = identifiers_mappings.get(db_ns, db_ns.lower())
    if identifiers_ns in identifiers_registry:
        if re.match(identifiers_registry[identifiers_ns]['pattern'], db_id):
            return
        else:
            raise InvalidIdentifier(f'{db_ns}:{db_id}')
    elif db_ns in non_registry or db_ns in non_grounding:
        return
    else:
        raise InvalidIdentifier(f'{db_ns}:{db_id}')


def validate_db_refs(db_refs):
    """Return True if all the entries in the given db_refs are valid.

    Parameters
    ----------
    db_refs : dict
        A dict of database references, typically part of an INDRA Agent.

    Returns
    -------
    bool
        True if all the entries are valid, else False.
    """
    try:
        assert_valid_db_refs(db_refs)
        return True
    except ValueError:
        return False


def assert_valid_db_refs(db_refs):
    """Raise InvalidIdentifier error if any of the entries in the given
    db_refs are invalid.

    Parameters
    ----------
    db_refs : dict
        A dict of database references, typically part of an INDRA Agent.
    """
    for db_ns, db_id in db_refs.items():
        assert_valid_ns(db_ns)
        assert_valid_id(db_ns, db_id)


def validate_statement(stmt):
    """Return True if all the groundings in the given statement are valid.

    Parameters
    ----------
    stmt : indra.statements.Statement
        An INDRA Statement to validate.

    Returns
    -------
    bool
        True if all the db_refs entries of the Agents in the given
        Statement are valid, else False.
    """
    try:
        assert_valid_statement(stmt)
        return True
    except ValueError:
        return False


def assert_valid_statement(stmt):
    """Raise InvalidIdentifier error if any of the groundings in the given
    statement are invalid.

    Parameters
    ----------
    stmt : indra.statements.Statement
        An INDRA Statement to validate.
    """
    for agent in stmt.real_agent_list():
        assert_valid_db_refs(agent.db_refs)


def assert_valid_text_refs(text_refs):
    """Return True if the given text refs are valid"""
    for ns in text_refs:
        if ns != ns.upper():
            raise InvalidTextRefs(ns)

    for ns, pattern in text_ref_patterns.items():
        if ns in text_refs:
            if not re.match(pattern, text_refs[ns]):
                raise InvalidTextRefs(f'{ns}:{text_refs[ns]}')


def assert_valid_pmid_text_refs(evidence):
    """Return True if the pmid attribute is consistent with text refs"""
    tr_pmid = evidence.text_refs.get('PMID')
    if tr_pmid is not None:
        if evidence.pmid != tr_pmid:
            raise InvalidTextRefs(f'Evidence pmid {evidence.pmid} doesn\'t '
                                  f'match text refs pmid {tr_pmid}')


def assert_valid_bio_context(context):
    # We shouldn't make a context if the attributes are all None
    if all(getattr(context, attr) is None for attr in BioContext.attrs):
        raise InvalidContext('All context attributes are None.')
    # Here we check db_refs validity
    for attr in BioContext.attrs:
        val = getattr(context, attr)
        if val is None:
            continue
        if val is not None and not isinstance(val, RefContext):
            raise InvalidContext(f'Invalid context entry for {attr}')
        assert_valid_db_refs(val.db_refs)


def assert_valid_context(context):
    if context is None:
        return
    elif isinstance(context, BioContext):
        assert_valid_bio_context(context)
    elif isinstance(context, WorldContext):
        return


def assert_valid_evidence(evidence):
    assert_valid_pmid_text_refs(evidence)
    assert_valid_text_refs(evidence.text_refs)
    assert_valid_context(evidence.context)


def validate_evidence(evidence):
    try:
        assert_valid_evidence(evidence)
        return True
    except ValueError:
        return False
