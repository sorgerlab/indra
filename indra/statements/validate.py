"""This module implements a number of functions that can be used to
validate INDRA Statements. The available functions include ones that
raise custom exceptions derived from ValueError if an invalidity is
found. These come with a helpful error message that can be caught
and printed to learn about the specific issue. Another set of functions
do not raise exceptions, rather, return True or False depending on whether
the given input is valid or invalid.

For validating namespaces and identifiers, there are two validators available,
one that uses data from identifiers.org and another for Bioregistry.
"""

import re
import logging
from indra.statements import *
from indra.databases import bioregistry, identifiers


logger = logging.getLogger(__name__)


# These are namespaces that can appear in db_refs but are actually not
# representing grounding.
non_grounding = {
    'TEXT', 'TEXT_NORM'
}


text_ref_patterns = {
    'PMID': re.compile(r'(\d+)'),
    'PMCID': re.compile(r'PMC(\d+)'),
    # See https://www.crossref.org/blog/dois-and-matching-regular-expressions/
    # here I added a-z to allow non-capitalized DOI text, we could
    # change this is we want strict capital letters in the DOI
    'DOI':  re.compile(r'^10.\d{4,9}/[-._;()/:A-Za-z0-9]+'),
}


class UnknownNamespace(ValueError):
    def __str__(self):
        return f'Unknown namespace: {self.args[0]}'


class MissingIdentifier(ValueError):
    """Raised when the identifier is None."""

    def __str__(self):
        return f'Missing identifier for {self.args[0]}'


class InvalidIdentifier(ValueError):
    """Raised when the identifier doesn't match the pattern."""

    def __str__(self):
        return f'Invalid identifier: {self.args[1]} for {self.args[0]} pattern {self.args[2]}'


class UnknownIdentifier(ValueError):
    """Raise when the database is neither registered with identifiers.org or
    manually added to the :data:`indra.databases.identifiers.non_registry` list.
    """

    def __str__(self):
        return f'Unknown identifier: {self.args[0]}:{self.args[1]}'


class InvalidTextRefs(ValueError):
    pass


class InvalidContext(ValueError):
    pass


class InvalidAgent(ValueError):
    pass


class InvalidStatement(ValueError):
    pass


class IdentifiersValidator:
    """A class that can be used to validate INDRA Statements."""

    @staticmethod
    def assert_valid_ns(db_ns):
        if db_ns in non_grounding or db_ns in identifiers.non_registry:
            return
        identifiers_ns = \
            identifiers.identifiers_mappings.get(db_ns, db_ns.lower())
        if identifiers_ns in identifiers.identifiers_registry:
            return
        raise UnknownNamespace(db_ns)

    @staticmethod
    def assert_valid_id(db_ns, db_id):
        if db_id is None:
            raise MissingIdentifier(db_ns, None)
        if db_ns in non_grounding or db_ns in identifiers.non_registry:
            return
        identifiers_ns = \
            identifiers.identifiers_mappings.get(db_ns, db_ns.lower())
        if identifiers_ns in identifiers.identifiers_registry:
            pattern = identifiers.identifiers_registry[
                identifiers_ns]['pattern_compiled']
            if pattern.match(db_id):
                return
            else:
                raise InvalidIdentifier(db_ns, db_id, pattern.pattern)
        else:
            raise UnknownIdentifier(db_ns, db_id)


class BioregistryValidator:
    """A class that can be used to validate INDRA Statements."""

    @staticmethod
    def assert_valid_ns(db_ns):
        if db_ns in non_grounding:
            return
        prefix = bioregistry.bioregistry_overrides.get(db_ns, db_ns.lower())
        if prefix in bioregistry.registry:
            return
        raise UnknownNamespace(db_ns)

    @staticmethod
    def assert_valid_id(db_ns, db_id):
        if db_id is None:
            raise MissingIdentifier(db_ns, None)
        if db_ns in non_grounding:
            return
        prefix = bioregistry.bioregistry_overrides.get(db_ns, db_ns.lower())
        if prefix in bioregistry.registry:
            pattern = bioregistry.registry[prefix].get('pattern_compiled')
            if not pattern or pattern.match(db_id):
                return
            else:
                raise InvalidIdentifier(db_ns, db_id, pattern.pattern)
        else:
            raise UnknownIdentifier(db_ns, db_id)


default_validator = IdentifiersValidator()


def validate_ns(db_ns, validator=default_validator):
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
        assert_valid_ns(db_ns, validator=validator)
        return True
    except ValueError:
        return False


def assert_valid_ns(db_ns, validator=default_validator):
    """Raise UnknownNamespace error if the given namespace is unknown.

    Parameters
    ----------
    db_ns : str
        The namespace.
    """
    return validator.assert_valid_ns(db_ns)


def validate_id(db_ns, db_id, validator=default_validator):
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
        assert_valid_id(db_ns, db_id, validator=validator)
        return True
    except ValueError:
        return False


def assert_valid_id(db_ns, db_id, validator=default_validator):
    """Raise InvalidIdentifier error if the ID is invalid in the given
    namespace.

    Parameters
    ----------
    db_ns : str
        The namespace.
    db_id : str
        The ID.
    """
    return validator.assert_valid_id(db_ns, db_id)


def validate_db_refs(db_refs, validator=default_validator):
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
        assert_valid_db_refs(db_refs, validator=validator)
        return True
    except ValueError:
        return False


def assert_valid_db_refs(db_refs, validator=default_validator):
    """Raise InvalidIdentifier error if any of the entries in the given
    db_refs are invalid.

    Parameters
    ----------
    db_refs : dict
        A dict of database references, typically part of an INDRA Agent.
    """
    for db_ns, db_id in db_refs.items():
        assert_valid_ns(db_ns, validator=validator)
        assert_valid_id(db_ns, db_id, validator=validator)


def assert_valid_agent(agent, validator=default_validator):
    """Raise InvalidAgent is there is an invalidity in the Agent.

    Parameters
    ----------
    agent : indra.statements.Agent
        The agent to check.
    """
    if agent is None:
        return
    if agent.name is None:
        raise InvalidAgent('Agent missing name')
    assert_valid_db_refs(agent.db_refs, validator=validator)


def validate_agent(agent, validator=default_validator):
    """Return False if is there is an invalidity in the Agent, otherwise True.

    Parameters
    ----------
    agent : indra.statements.Agent
        The agent to check.

    Returns
    -------
    bool
        True if the agent is valid, False otherwise.
    """
    try:
        assert_valid_agent(agent, validator=validator)
        return True
    except ValueError:
        return False


def assert_valid_statement_semantics(stmt):
    """Raise InvalidStatement error if the given statement is invalid.

    Parameters
    ----------
    statement : indra.statements.Statement
        The statement to check.
    """
    if all(a is None for a in stmt.agent_list()):
        raise InvalidStatement('Statement with all None agents')

    if isinstance(stmt, Complex):
        if any(m is None for m in stmt.members):
            raise InvalidStatement('Complex with None agent.')
    elif isinstance(stmt, Conversion):
        if any(o is None for o in stmt.obj_from):
            raise InvalidStatement('Conversion with None reactant.')
        if any(o is None for o in stmt.obj_to):
            raise InvalidStatement('Conversion with None product.')
    elif isinstance(stmt, RegulateActivity):
        if stmt.subj is None:
            raise InvalidStatement('Regulation missing subject.')
        if stmt.obj is None:
            raise InvalidStatement('Regulation missing object.')
    elif isinstance(stmt, RegulateAmount):
        if stmt.obj is None:
            raise InvalidStatement('Regulation missing object.')
    elif isinstance(stmt, (Gap, Gef)):
        if any(a is None for a in stmt.agent_list()):
            raise InvalidStatement('Gap/Gef missing agent.')
    elif isinstance(stmt, Modification):
        if stmt.sub is None:
            raise InvalidStatement('Modification missing substrate.')
    elif isinstance(stmt, Translocation):
        if stmt.from_location is None and stmt.to_location is None:
            raise InvalidStatement('Translocation with no locations.')


def validate_statement(stmt, validator=default_validator):
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
        assert_valid_statement(stmt, validator=validator)
        return True
    # Some deeper validation checks in statements raise TypeErrors
    except (ValueError, TypeError):
        return False


def assert_valid_statement(stmt, validator=default_validator):
    """Raise an error if there is anything invalid in the given statement.

    Parameters
    ----------
    stmt : indra.statements.Statement
        An INDRA Statement to validate.
    """
    assert_valid_statement_semantics(stmt)
    for agent in stmt.real_agent_list():
        assert_valid_agent(agent, validator=validator)
    for ev in stmt.evidence:
        assert_valid_evidence(ev, validator=validator)


def assert_valid_statements(stmts, validator=default_validator):
    """Raise an error of any of the given statements is invalid.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements to validate.
    """
    for stmt in stmts:
        assert_valid_statement(stmt, validator=validator)


def print_validation_report(stmts, validator=default_validator):
    """Log the first validation error encountered for each given statement.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements to validate.
    """
    for idx, stmt in enumerate(stmts):
        try:
            assert_valid_statement(stmt, validator=validator)
        except Exception as e:
            logger.info(f'{idx}: {type(e).__name__} - {e}')


def assert_valid_text_refs(text_refs):
    """Raise an InvalidTextRefs error if the given text refs are invalid."""
    for ns in text_refs:
        if ns != ns.upper():
            raise InvalidTextRefs(ns)

    for ns, pattern in text_ref_patterns.items():
        if ns in text_refs:
            if text_refs[ns] is None:
                raise InvalidTextRefs(f'{ns}:{text_refs[ns]}')
            elif not re.match(pattern, text_refs[ns]):
                raise InvalidTextRefs(f'{ns}:{text_refs[ns]}')


def validate_text_refs(text_refs):
    """Return True if the given text refs are valid."""
    try:
        assert_valid_text_refs(text_refs)
        return True
    except ValueError:
        return False


def assert_valid_pmid_text_refs(evidence):
    """Return True if the pmid attribute is consistent with text refs"""
    tr_pmid = evidence.text_refs.get('PMID')
    if tr_pmid is not None:
        if evidence.pmid != tr_pmid:
            raise InvalidTextRefs(f'Evidence pmid {evidence.pmid} doesn\'t '
                                  f'match text refs pmid {tr_pmid}')


def assert_valid_bio_context(context, validator=default_validator):
    """Raise InvalidContext error if the given bio-context is invalid.

    Parameters
    ----------
    context : indra.statements.BioContext
        The context object to validate.
    """
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
        assert_valid_db_refs(val.db_refs, validator=validator)


def assert_valid_context(context, validator=default_validator):
    """Raise InvalidContext error if the given context is invalid.

    Parameters
    ----------
    context : indra.statements.Context
        The context object to validate.
    """
    if context is None:
        return
    elif isinstance(context, BioContext):
        assert_valid_bio_context(context, validator=validator)
    elif isinstance(context, WorldContext):
        return


def assert_valid_evidence(evidence, validator=default_validator):
    """Raise an error if the given evidence is invalid.

    Parameters
    ----------
    evidence : indra.statements.Evidence
        The evidence object to validate.
    """
    assert_valid_pmid_text_refs(evidence)
    assert_valid_text_refs(evidence.text_refs)
    assert_valid_context(evidence.context, validator=validator)


def validate_evidence(evidence, validator=default_validator):
    """Return False if the given evidence is invalid, otherwise True.

    Parameters
    ----------
    evidence : indra.statements.Evidence
        The evidence object to validate.

    Returns
    -------
    bool
        True if the evidence is valid, otherwise False.
    """
    try:
        assert_valid_evidence(evidence, validator=validator)
        return True
    except ValueError:
        return False
