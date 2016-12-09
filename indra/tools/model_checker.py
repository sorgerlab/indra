from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging
import numpy as np
import networkx
import itertools
from copy import deepcopy
from collections import deque
from itertools import product
from pysb import kappa
from pysb import Observable, ComponentSet
from pysb.core import as_complex_pattern
from indra.statements import *
from indra.assemblers import pysb_assembler as pa
logger = logging.getLogger('model_checker')

class ModelChecker(object):
    """Check a PySB model against a set of INDRA statements."""

    def __init__(self, model, stmts_to_check=None):
        self.model = model
        if stmts_to_check:
            self.statements = stmts_to_check
        else:
            self.statements = []
        self._im = None

    def get_im(self, force_update=True):
        if self._im and not force_update:
            return self._im
        if not self.model:
            raise Exception("Cannot get influence map if there is no model.")
        else:
            self._im = kappa.influence_map(self.model)
            self._im.is_multigraph = lambda: False
            return self._im

    def check_model(self):
        results = []
        for stmt in self.statements:
            result = self.check_statement(stmt)
            results.append((stmt, result))
        return results

    def check_statement(self, stmt):
        if isinstance(stmt, Modification):
            return self.check_modification(stmt)
        elif isinstance(stmt, Activation):
            return self.check_activation(stmt)
        else:
            return False

    def check_activation(self, stmt):
        logger.info('Checking stmt: %s' % stmt)
        # FIXME Currently this will match rules with the corresponding monomer
        # pattern from the Activation statement, which will nearly always have
        # no state conditions on it. In future, this statement should also match
        # rules in which 1) the agent is in its active form, or 2) the agent is
        # tagged as the enzyme in a rule of the appropriate activity (e.g., a
        # phosphorylation rule) FIXME
        subj_mp = pa.get_monomer_pattern(self.model, stmt.subj)
        target_polarity = 1 if stmt.is_activation else -1
        # This may fail, since there may be no rule in the model activating the
        # object, and the object may not have an "active" site of the
        # appropriate type
        obj_obs_name = pa.get_agent_rule_str(stmt.obj) + '_obs'
        try:
            obj_site_pattern = pa.get_site_pattern(stmt.obj)
            obj_site_pattern.update({pa.active_site_names[stmt.obj_activity]:
                                     'active'})
            obj_monomer = self.model.monomers[stmt.obj.name]
            obj_mp = obj_monomer(**obj_site_pattern)
        except Exception as e:
            logger.info("Could not create obj monomer pattern: %s" % e)
            return False
        obj_obs = Observable(obj_obs_name, obj_mp, _export=False)
        return self._find_im_paths(subj_mp, obj_obs, target_polarity)

    # Need to be able to take an arbitrary set of agent modification/mut etc.
    # conditions and generate an appropriate set of monomer patterns and/or
    # observables. This has to apply to both input rules/agents as well as
    # targets. Imagine the case of
    # Phosphorylation(MAP2K1(mods:[phosphorylation, None, None]), MAPK1)
    # Here, we need to match any rule where MEK (in any phosphorylated state)
    # phosphoryaltes ERK on any site.
    # If this can be done, then resolving activations becomes simply a matter of
    # looking up which states are active, just like looking up phosphorylations.
    # So we'll start by looking through the modification conditions of an agent.

    def check_modification(self, stmt):
        # Identify the observable we're looking for in the model, which
        # may not exist!
        # The observable is the modified form of the substrate
        logger.info('Checking stmt: %s' % stmt)
        # Look for an agent with the appropriate grounding in the model
        if stmt.enz is not None:
            enz_mps = list(pa.grounded_monomer_patterns(self.model, stmt.enz))
            if not enz_mps:
                logger.info('No monomers found corresponding to agent %s' %
                             stmt.enz)
                return False
        else:
            enz_mps = [None]
        # Get target polarity
        demodify_list = (Dephosphorylation, Dehydroxylation, Desumoylation,
                         Deacetylation, Deglycosylation, Deribosylation,
                         Deubiquitination, Defarnesylation)
        target_polarity = -1 if type(stmt) in demodify_list else 1
        # Add modification to substrate agent
        modified_sub = _add_modification_to_agent(stmt.sub, 'phosphorylation',
                                                  stmt.residue, stmt.position)
        # Now look up the modified agent in the model
        #obj_mp = pa.get_monomer_pattern(self.model, modified_sub)
        #site_pattern = pa.get_site_pattern(modified_sub)
        #monomer = self.model.monomers[modified_sub.name]
        #try:
        #    obj_mp = monomer(**site_pattern)
        #except Exception as e:
        #    logger.info("Invalid site: %s" % e)
        #    return False
        obs_name = pa.get_agent_rule_str(modified_sub) + '_obs'
        obj_mps = list(pa.grounded_monomer_patterns(self.model, modified_sub))
        if not obj_mps:
            logger.info('Failed to create observable; returning False')
            return False
        for enz_mp, obj_mp in itertools.product(enz_mps, obj_mps):
            obj_obs = Observable(obs_name, obj_mp, _export=False)
            if self._find_im_paths(enz_mp, obj_obs, target_polarity):
                return True
        # If we got here, then there was no path for any observable
        return False

    def _find_im_paths(self, subj_mp, obj_obs, target_polarity):
        # Reset the observables list
        self.model.observables = ComponentSet([])
        self.model.add_component(obj_obs)
        # Find rules in the model corresponding to the input
        logger.info('Finding paths between %s and %s with polarity %s' %
                    (subj_mp, obj_obs, target_polarity))
        if subj_mp is None:
            input_rule_set = None
        else:
            input_rules = match_lhs(subj_mp, self.model.rules)
            input_rule_set = set([r.name for r in input_rules])
                                  #if r.name.startswith(stmt.enz.name)])
            logger.info('Found %s input rules matching %s' %
                        (len(input_rules), str(subj_mp)))
            # If we have enzyme information but there are no input rules
            # matching the enzyme, then there is no path
            if not input_rules:
                return False
        # Generate the predecessors to our observable
        num_paths = 0
        for (source, polarity, path_length) in \
                    _find_sources(self.get_im(), obj_obs.name, input_rule_set,
                                  target_polarity):
            num_paths += 1
        if num_paths > 0:
            return True
        else:
            return False


def _add_modification_to_agent(agent, mod_type, residue, position):
    new_mod = ModCondition(mod_type, residue, position)
    # Check if this modification already exists
    for old_mod in agent.mods:
        if old_mod.equals(new_mod):
            return agent
    new_agent = deepcopy(agent)
    new_agent.mods.append(new_mod)
    return new_agent


def match_lhs(cp, rules):
    rule_matches = []
    for rule in rules:
        reactant_pattern = rule.rule_expression.reactant_pattern
        for rule_cp in reactant_pattern.complex_patterns:
            if cp_embeds_into(rule_cp, cp):
                rule_matches.append(rule)
                break
    return rule_matches


def match_rhs(cp, rules):
    rule_matches = []
    for rule in rules:
        product_pattern = rule.rule_expression.product_pattern
        for rule_cp in product_pattern.complex_patterns:
            if cp_embeds_into(rule_cp, cp):
                rule_matches.append(rule)
                break
    return rule_matches

"""
def find_production_rules(cp, rules):
    # Find rules where the CP matches the left hand side
    lhs_rule_set = set(match_lhs(cp, rules))
    # Now find rules where the CP matches the right hand side
    rhs_rule_set = set(match_rhs(cp, rules))
    # Production rules are rules where there is a match on the right hand
    # side but not on the left hand side
    prod_rules = list(rhs_rule_set.difference(lhs_rule_set))
    return prod_rules

def find_consumption_rules(cp, rules):
    # Find rules where the CP matches the left hand side
    lhs_rule_set = set(match_lhs(cp, rules))
    # Now find rules where the CP matches the right hand side
    rhs_rule_set = set(match_rhs(cp, rules))
    # Consumption rules are rules where there is a match on the left hand
    # side but not on the right hand side
    cons_rules = list(lhs_rule_set.difference(rhs_rule_set))
    return cons_rules
"""


def cp_embeds_into(cp1, cp2):
    # Check that any state in cp2 is matched in cp2
    # If the thing we're matching to is just a monomer pattern, that makes
    # things easier--we just need to find the corresponding monomer pattern
    # in cp1
    cp1 = as_complex_pattern(cp1)
    cp2 = as_complex_pattern(cp2)
    if len(cp2.monomer_patterns) == 1:
        mp2 = cp2.monomer_patterns[0]
        # Iterate over the monomer patterns in cp1 and see if there is one
        # that has the same name
        for mp1 in cp1.monomer_patterns:
            if mp_embeds_into(mp1, mp2):
                return True
    return False


def mp_embeds_into(mp1, mp2):
    sc_matches = []
    if mp1.monomer.name != mp2.monomer.name:
        return False
    # Check that all conditions in mp2 are met in mp1
    for site_name, site_state in mp2.site_conditions.items():
        if site_name not in mp1.site_conditions or \
           site_state != mp1.site_conditions[site_name]:
            return False
    return True


def _get_signed_predecessors(im, node, polarity):
    signed_pred_list = []
    predecessors = im.predecessors_iter
    for pred in predecessors(node):
        pred_edge = im.get_edge(pred, node)
        yield (pred, _get_edge_sign(pred_edge) * polarity)


def _find_sources(im, target, sources, polarity):
    # First, generate a list of visited nodes
    # Copied/adapted from networkx
    visited = set([(target, 1)])
    # Generate list of predecessor nodes with a sign updated according to the
    # sign of the target node
    target_tuple = (target, 1)
    # The queue holds tuples of "parents" (in this case downstream nodes) and
    # their "children" (in this case their upstream influencers)
    queue = deque([(target_tuple, _get_signed_predecessors(im, target, 1), 1)])
    while queue:
        parent, children, path_length = queue[0]
        try:
            # Get the next child in the list
            (child, sign) = next(children)
            # The sign of the child should be updated according to the sign of
            # the parent!
            # Could do this in the expansion step itself, or here;
            # If done during the expansion step, 
            # Check this child against the visited list. If we haven't visited
            # it, then we return the parent/child pair
            if (sources is None or child in sources) and sign == polarity:
                logger.info("Found path to %s from %s with desired sign %s "
                            "with length %d" %
                            (target, child, polarity, path_length))
                yield (child, sign, path_length)
            if (child, sign) not in visited:
                # TODO: Update this to include the path length, and perhaps to
                # check for sources matching the list
                visited.add((child, sign))
                queue.append((child, _get_signed_predecessors(im, child, sign),
                              path_length + 1))
                # Add in the new parent child pair after expanding out
        # Once we've finished iterating over the children of the current node,
        # pop the node off and go to the next one in the queue
        except StopIteration:
            queue.popleft()
            path_length += 1
    return None


def positive_path(im, path):
    # This doesn't address the effect of the rules themselves on the
    # observables of interest--just the effects of the rules on each other
    edge_polarities = []
    for rule_ix in range(len(path) - 1):
        from_rule = path[rule_ix]
        to_rule = path[rule_ix + 1]
        edge = im.get_edge(from_rule, to_rule)
        edge_polarities.append(_get_edge_sign(edge))
    # Compute and return the overall path polarity
    path_polarity = np.prod(edge_polarities)
    assert path_polarity == 1 or path_polarity == -1
    return True if path_polarity == 1 else False


def _get_edge_sign(edge):
    if edge.attr.get('color') is None:
        raise Exception('No color attribute for edge.')
    elif edge.attr['color'] == 'green':
        return 1
    elif edge.attr['color'] == 'red':
        return -1
    else:
        raise Exception('Unexpected edge color: %s' % edge.attr['color'])


