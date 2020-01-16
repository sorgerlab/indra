"""This module implements a client to the Gilda grounding web service,
and contains functions to help apply it during the course of INDRA assembly."""
import os
import logging
import requests
from urllib.parse import urljoin
from indra.preassembler.grounding_mapper.standardize \
    import standardize_agent_name
from indra.config import get_config, has_config
from .adeft import _get_text_for_grounding

logger = logging.getLogger(__name__)

grounding_service_url = get_config('GILDA_URL') if has_config('GILDA_URL') \
    else 'http://grounding.indra.bio'


def get_gilda_models(offline=True):
    """Return a list of strings for which Gilda has a disambiguation model.

    Parameters
    ----------
    offline : Optional[bool]
        If True, INDRA's local resource file is used to determine the set
        of models. Otherwise, the web service is invoked to get the
        set of models. Default: True

    Returns
    -------
    list[str]
        A list of entity strings.
    """
    if offline:
        fname = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                             os.pardir, os.pardir, 'resources',
                             'gilda_models.txt')
        with open(fname, 'r') as fh:
            return [l.strip() for l in fh.readlines()]
    else:
        res = requests.post(urljoin(grounding_service_url, 'models'))
        models = res.json()
        return models


def ground_agent(agent, txt, context=None):
    resp = requests.post(urljoin(grounding_service_url, 'ground'),
                         json={'text': txt, 'context': context})
    results = resp.json()
    if results:
        db_refs = {'TEXT': txt,
                   results[0]['term']['db']: results[0]['term']['id']}
        agent.db_refs = db_refs
        standardize_agent_name(agent, standardize_refs=True)
    return results


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
            ground_agent(agent, txt, context)


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


def run_gilda_disambiguation(stmt, agent, idx):
    """Run Gilda disambiguation on an Agent in a given Statement.

    This function looks at the evidence of the given Statement and attempts
    to look up the full paper or the abstract for the evidence. If both of
    those fail, the evidence sentence itself is used for disambiguation.
    The disambiguation model corresponding to the Agent text is then called,
    and the highest scoring returned grounding is set as the Agent's new
    grounding.

    The Statement's annotations as well as the Agent are modified in place
    and no value is returned.

    Parameters
    ----------
    stmt : indra.statements.Statement
        An INDRA Statement in which the Agent to be disambiguated appears.
    agent : indra.statements.Agent
        The Agent (potentially grounding mapped) which we want to
        disambiguate in the context of the evidence of the given Statement.
    idx : int
        The index of the new Agent's position in the Statement's agent list
        (needed to set annotations correctly).

    Returns
    -------
    bool
        True if disambiguation was successfully applied, and False otherwise.
        Reasons for a False response can be the lack of evidence as well as
        failure to obtain text for grounding disambiguation.
    """
    success = False
    # If the Statement doesn't have evidence for some reason, then there is
    # no text to disambiguate by
    # NOTE: we might want to try disambiguating by other agents in the
    # Statement
    if not stmt.evidence:
        return False
    # Initialize annotations if needed so predicted
    # probabilities can be added to Agent annotations
    annots = stmt.evidence[0].annotations
    agent_txt = agent.db_refs['TEXT']
    if 'agents' in annots:
        if 'gilda' not in annots['agents']:
            annots['agents']['gilda'] = \
                {'gilda': [None for _ in stmt.agent_list()]}
    else:
        annots['agents'] = {'gilda': [None for _ in stmt.agent_list()]}
    grounding_text = _get_text_for_grounding(stmt, agent_txt)
    if grounding_text:
        gilda_result = ground_agent(agent, agent_txt, grounding_text)
        if gilda_result:
            annots['agents']['gilda'][idx] = gilda_result
            success = True
    return success
