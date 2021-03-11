"""This module implements a client to the Gilda grounding web service,
and contains functions to help apply it during the course of INDRA assembly."""

import logging
import requests
from copy import deepcopy
from typing import Any, Callable, List, Mapping, Optional, Tuple
from urllib.parse import urljoin
from indra.ontology.standardize \
    import standardize_agent_name
from indra.config import get_config, has_config
from indra.pipeline import register_pipeline


logger = logging.getLogger(__name__)

grounding_service_url = get_config('GILDA_URL', failure_ok=True) \
    if has_config('GILDA_URL') else 'http://grounding.indra.bio/'


def get_grounding(
    txt: str,
    context: Optional[str] = None,
    mode: Optional[str] = 'web',
) -> Tuple[Mapping[str, Any], List[Any]]:
    """Return the top Gilda grounding for a given text.

    Parameters
    ----------
    txt : str
        The text to ground.
    context : Optional[str]
        Any context for the grounding.
    mode : Optional[str]
        If 'web', the web service given in the GILDA_URL config setting or
        environmental variable is used. Otherwise, the gilda package is
        attempted to be imported and used. Default: web

    Returns
    -------
    dict
        If no grounding was found, it is an empty dict. Otherwise, it's a
        dict with the top grounding returned from Gilda.
    list
        The list of ScoredMatches
    """
    grounding = {}
    if mode == 'web':
        resp = requests.post(urljoin(grounding_service_url, 'ground'),
                             json={'text': txt, 'context': context})
        results = resp.json()
        if results:
            grounding = {results[0]['term']['db']: results[0]['term']['id']}
    else:
        from gilda import ground
        results = ground(txt, context)
        if results:
            grounding = {results[0].term.db: results[0].term.id}
            results = [sm.to_json() for sm in results]
    return grounding, results


def get_gilda_models(mode='web'):
    """Return a list of strings for which Gilda has a disambiguation model.

    Parameters
    ----------
    mode : Optional[str]
        If 'web', the web service given in the GILDA_URL config setting or
        environmental variable is used. Otherwise, the gilda package is
        attempted to be imported and used. Default: web

    Returns
    -------
    list[str]
        A list of entity strings.
    """
    if mode == 'web':
        res = requests.post(urljoin(grounding_service_url, 'models'))
        models = res.json()
        return models
    else:
        from gilda import get_models
        return get_models()


def ground_agent(agent, txt, context=None, mode='web'):
    """Set the grounding of a given agent, by re-grounding with Gilda.

    This function changes the agent in place without returning a value.

    Parameters
    ----------
    agent : indra.statements.Agent
        The Agent whose db_refs shuld be changed.
    txt : str
        The text by which the Agent should be grounded.
    context : Optional[str]
        Any additional text context to help disambiguate the sense
        associated with txt.
    mode : Optional[str]
        If 'web', the web service given in the GILDA_URL config setting or
        environmental variable is used. Otherwise, the gilda package is
        attempted to be imported and used. Default: web
    """
    gr, results = get_grounding(txt, context, mode)
    if gr:
        db_refs = {'TEXT': txt}
        db_refs.update(gr)
        agent.db_refs = db_refs
        standardize_agent_name(agent, standardize_refs=True)
    return results


def ground_statement(stmt, mode='web', ungrounded_only=False):
    """Set grounding for Agents in a given Statement using Gilda.

    This function modifies the original Statement/Agents in place.

    Parameters
    ----------
    stmt : indra.statements.Statement
        A Statement to ground
    mode : Optional[str]
        If 'web', the web service given in the GILDA_URL config setting or
        environmental variable is used. Otherwise, the gilda package is
        attempted to be imported and used. Default: web
    ungrounded_only : Optional[str]
        If True, only ungrounded Agents will be grounded, and ones that
        are already grounded will not be modified. Default: False
    """
    if stmt.evidence and stmt.evidence[0].text:
        context = stmt.evidence[0].text
    else:
        context = None
    for agent in stmt.agent_list():
        if agent is not None and 'TEXT' in agent.db_refs:
            txt = agent.db_refs['TEXT']
            gr = agent.get_grounding()
            if not ungrounded_only or gr[0] is None:
                ground_agent(agent, txt, context, mode=mode)


@register_pipeline
def ground_statements(stmts, mode='web', sources=None, ungrounded_only=False):
    """Set grounding for Agents in a list of Statements using Gilda.

    This function modifies the original Statements/Agents in place.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of Statements to ground
    mode : Optional[str]
        If 'web', the web service given in the GILDA_URL config setting or
        environmental variable is used. Otherwise, the gilda package is
        attempted to be imported and used. Default: web
    sources : Optional[list]
        If given, only statements from the given sources are grounded. The
        sources have to correspond to valid source_api entries, e.g.,
        'reach', 'sparser', etc. If not given, statements from all
        sources are grounded.
    ungrounded_only : Optional[str]
        If True, only ungrounded Agents will be grounded, and ones that
        are already grounded will not be modified. Default: False

    Returns
    -------
    list[indra.statement.Statements]
        The list of Statements that were changed in place by reference.
    """
    source_filter = set(sources) if sources else set()
    grounded_stmts = deepcopy(stmts)
    for stmt in grounded_stmts:
        if not source_filter or (stmt.evidence and stmt.evidence[0].source_api
                                 in source_filter):
            ground_statement(stmt, mode=mode, ungrounded_only=ungrounded_only)
    return grounded_stmts
