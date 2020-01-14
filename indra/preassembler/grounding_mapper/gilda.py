"""This module implements a client to the Gilda grounding web service,
and contains functions to help apply it during the course of INDRA assembly."""

import requests
from .mapper import GroundingMapper

grounding_service_url = 'http://grounding.indra.bio'


def get_gilda_models():
    """Return a list of strings for which Gilda has a disambiguation model.

    Returns
    -------
    list[str]
        A list of entity strings.
    """
    res = requests.post(grounding_service_url + '/models')
    models = res.json()
    return models


def ground_statement(stmt):
    """Set grounding for Agents in a given Statement using Gilda.

    This function modifies the original Statement/Agents in place.

    Parameters
    ----------
    stmt : indra.statements.Statement
        A Statement to ground
    """
    if stmt.evidence and stmt.evidence[0].text:
        context = stmt.evidence[0].text
    else:
        context = None
    for agent in stmt.agent_list():
        if agent is not None and 'TEXT' in agent.db_refs:
            txt = agent.db_refs['TEXT']
            resp = requests.post(grounding_service_url + '/ground',
                                 json={'text': txt,
                                       'context': context})
            results = resp.json()
            if results:
                db_refs = {'TEXT': txt,
                           results[0]['term']['db']:
                               results[0]['term']['id']}
                agent.db_refs = db_refs
                GroundingMapper.standardize_agent_name(agent,
                                                       standardize_refs=True)


def ground_statements(stmts):
    """Set grounding for Agents in a list of Statements using Gilda.

    This function modifies the original Statements/Agents in place.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of Statements to ground
    """
    for stmt in stmts:
        ground_statement(stmt)
