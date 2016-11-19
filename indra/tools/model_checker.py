from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging
import numpy as np
import networkx
import itertools
from copy import deepcopy
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
        else:
            return False

    def check_modification(self, stmt):
        # Identify the observable we're looking for in the model, which
        # may not exist!
        # The observable is the modified form of the substrate
        logger.info('Checking stmt: %s' % stmt)
        # Look for an agent with the appropriate grounding in the model
        #enz_db_refs = stmt.enz.db_refs
        if stmt.enz is not None:
            enz_mp = pa.get_monomer_pattern(self.model, stmt.enz)
        else:
            enz_mp = None

        if type(stmt) in (Dephosphorylation,):
            unmodify_stmt = True
        else:
            unmodify_stmt = False
        modified_sub = _add_modification_to_agent(stmt.sub, 'phosphorylation',
                                                  stmt.residue, stmt.position)
        site_pattern = pa.get_site_pattern(modified_sub)
        obs_name = pa.get_agent_rule_str(modified_sub) + '_obs'
        monomer = self.model.monomers[modified_sub.name]
        # If we try to get a monomer pattern for a monomer that doesn't have
        # the specified site, PySB throws a (generic) Exception; we catch this
        # as an indication that the site doesn't exist, and hence the
        # modification is not observed, and return False.
        try:
            obs = Observable(obs_name, monomer(**site_pattern), _export=False)
        except Exception as e:
            logger.info("Invalid site: %s" % e)
            return False
        # Reset the observables list
        self.model.observables = ComponentSet([])
        self.model.add_component(obs)
        # Find rules in the model corresponding to the input
        if enz_mp is not None:
            input_rules = match_lhs(enz_mp, self.model.rules)
            input_rule_set = set([r.name for r in input_rules])
                                  #if r.name.startswith(stmt.enz.name)])
            logger.info('Found %s input rules matching %s' %
                        (len(input_rules), str(enz_mp)))
            # If we have enzyme information but there are no input rules matching
            # the enzyme, then there is no path
            if not input_rules:
                return False
        # Generate the predecessors to our observable
        pred_dict = networkx.bfs_predecessors(self.get_im(), obs_name)
        preds = set(pred_dict.keys())
        # If we're missing participant A (no enzyme), then we only need to check
        # that there are any predecessors to our observable for this statement
        # to be valid
        if enz_mp is None and preds:
            return True
        elif enz_mp is None:
            return False
        # If we've got this far, then we know that we have an enzyme (enz_mp is
        # not None), and there are input rules matching this enzyme
        # (input_rule_set is not empty)
        assert enz_mp and input_rule_set
        # Check to see if any of our input rules are in our set of upstream
        # influencer rules
        source_influences = input_rule_set.intersection(set(preds))
        logger.info('Source influences')
        logger.info(source_influences)
        # If the set is non-empty, then we have a causal path
        if source_influences:
            return True
        # Otherwise, there is no path
        else:
            return False
        #logger.info("Trying input rule %s" % input_rule.name)
        #sp_gen = networkx.all_simple_paths(self.get_im(),
        #                                   input_rule.name,
        #                                   obs_name)
        #if networkx.has_path(self.get_im(), input_rule.name,
        #                     obs_name):
        #    logger.info('Found path between %s and %s' %
        #                (input_rule.name, obs_name))
        #    return True
        #else:
        #    logger.info('No path between %s and %s' %
        #                (input_rule.name, obs_name))
        #    continue
        # Iterate over paths until we find one with positive
        # polarity
        #need_positive_path = False if unmodify_stmt else True
        #for sp in sp_gen:
        #    logger.info('Found path: %s' % str(sp))
        #    if need_positive_path == \
        #            positive_path(self.get_im(), sp):
        #        logger.info('Found non-repeating path, '
        #                    'length %d: %s' %
        #                    (len(sp), str(sp)))
        #        paths.append(sp)
        #        break


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


def positive_path(im, path):
    # This doesn't address the effect of the rules themselves on the
    # observables of interest--just the effects of the rules on each other
    edge_polarities = []
    for rule_ix in range(len(path) - 1):
        from_rule = path[rule_ix]
        to_rule = path[rule_ix + 1]
        edge = im.get_edge(from_rule, to_rule)
        if _is_positive_edge(edge):
            edge_polarities.append(1)
        else:
            edge_polarities.append(-1)
    # Compute and return the overall path polarity
    path_polarity = np.prod(edge_polarities)
    assert path_polarity == 1 or path_polarity == -1
    return True if path_polarity == 1 else False


def _is_positive_edge(edge):
    if edge.attr.get('color') is None:
        raise Exception('No color attribute for edge.')
    elif edge.attr['color'] == 'green':
        return True
    elif edge.attr['color'] == 'red':
        return False
    else:
        raise Exception('Unexpected edge color: %s' % edge.attr['color'])


