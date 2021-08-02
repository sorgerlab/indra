"""This module contains a library of functions that are
useful for building custom preassembly logic for some applications.
They are typically used as matches_fun or refinement_fun arguments
to the Preassembler and other modules."""
from indra.statements import *
from indra.pipeline import register_pipeline


@register_pipeline
def agent_grounding_matches(agent):
    """Return an Agent matches key just based on grounding, not state."""
    if agent is None:
        return None
    return str(agent.entity_matches_key())


@register_pipeline
def agents_stmt_type_matches(stmt):
    """Return a matches key just based on Agent grounding and Stmt type."""
    agents = [agent_grounding_matches(a) for a in stmt.agent_list()]
    key = str((stmt.__class__.__name__, agents))
    return key


@register_pipeline
def agent_name_matches(agent):
    """Return a sorted, normalized bag of words as the name."""
    if agent is None:
        return None
    bw = '_'.join(sorted(list(set(agent.name.lower().split()))))
    return bw


@register_pipeline
def agent_name_stmt_type_matches(stmt):
    """Return True if the statement type and normalized agent name matches."""
    agents = [agent_name_matches(a) for a in stmt.agent_list()]
    key = str((stmt.__class__.__name__, agents))
    return key


@register_pipeline
def agent_name_stmt_matches(stmt):
    """Return the normalized agent names."""
    agents = [agent_name_matches(a) for a in stmt.agent_list()]
    key = str(agents)
    return key


@register_pipeline
def agent_name_polarity_matches(stmt):
    """Return a key for normalized agent names and polarity."""
    default_sign_dict = {'Activation': 0,
                         'Inhibition': 1,
                         'IncreaseAmount': 0,
                         'DecreaseAmount': 1}
    agents = [agent_name_matches(a) for a in stmt.agent_list()]
    pol = default_sign_dict[type(stmt).__name__]
    key = str((agents, pol))
    return key
