from pysb.core import as_complex_pattern
from pysb import kappa
from pysb import Observable
from indra.statements import *
from indra.assemblers import pysb_assembler as pa
from copy import deepcopy
import networkx
import itertools
import logging
import numpy as np

logger = logging.getLogger('model_checker')

class ModelChecker(object):
    """Check a PySB model against a set of INDRA statements."""

    def __init__(self, model, statements):
        self.model = model
        self.statements = statements
        self._im = None

    def get_im(self):
        if self._im:
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
        if isinstance(stmt, Phosphorylation):
            return self.check_phosphorylation(stmt)
        else:
            return False

    def check_phosphorylation(self, stmt):
        # Identify the observable we're looking for in the model, which
        # may not exist!
        # The observable is the modified form of the substrate
        enz_mp = pa.get_monomer_pattern(self.model, stmt.enz)
        modified_sub = _add_modification_to_agent(stmt.sub, 'phosphorylation',
                                                  stmt.residue, stmt.position)
        sub_mp = pa.get_monomer_pattern(self.model, modified_sub)
        # Generate the influence map
        # Find rules in the model corresponding to the inputs and outputs
        input_rules = match_lhs(enz_mp, self.model.rules)
        output_rules = match_rhs(sub_mp, self.model.rules)
        input_rule_names = [r.name for r in input_rules]
        output_rule_names = [r.name for r in output_rules]
        print input_rules
        print output_rules
        if input_rule_names and output_rule_names:
            # Generate the influence map
            for input_rule, output_rule in itertools.product(input_rule_names,
                                                             output_rule_names):
                try:
                    logger.info('Looking for path between %s and %s' %
                                (input_rule, output_rule))
                    sp_gen = networkx.shortest_simple_paths(self.get_im(),
                                                       input_rule,
                                                       output_rule)
                    # Iterate over paths until we find one with the desired
                    # polarity
                    for path_ix, sp in enumerate(sp_gen):
                        logger.info('Found simple path %s: %s' %
                                    (path_ix, str(sp)))
                        pol = path_polarity(self.get_im(), sp)
                        logger.info('Path %s has polarity %s' %
                                    (str(sp), pol))
                except networkx.NetworkXNoPath as nopath:
                    logger.info('No path found between %s and %s' %
                                (input_rule, output_rule))
            return False
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
    return rule_matches

def match_rhs(cp, rules):
    rule_matches = []
    for rule in rules:
        product_pattern = rule.rule_expression.product_pattern
        for rule_cp in product_pattern.complex_patterns:
            if cp_embeds_into(rule_cp, cp):
                rule_matches.append(rule)
    return rule_matches

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
    if len(path) == 1:
        return True
    edge_polarities = []
    for rule_ix in range(len(path) - 1):
        from_rule = path[rule_ix]
        to_rule = path[rule_ix + 1]
        edge = im.get_edge(from_rule, to_rule)
        if edge.attr.get('color') is None:
            raise Exception('No color attribute for edge.')
        elif edge.attr['color'] == 'green':
            edge_polarities.append(1)
        elif edge.attr['color'] == 'red':
            edge_polarities.append(-1)
        else:
            raise Exception('Unexpected edge color: %s' % edge.attr['color'])
    # Compute and return the overall path polarity
    path_polarity = np.prod(edge_polarities)
    assert path_polarity == 1 or path_polarity == -1
    return True if path_polarity == 1 else False
