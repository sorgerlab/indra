from pysb.core import as_complex_pattern
from pysb import kappa
from pysb import Observable
from indra.statements import *
from indra.assemblers import pysb_assembler as pa
from copy import deepcopy

class ModelChecker(object):
    """Check a PySB model against a set of INDRA statements."""

    def __init__(self, model, statements):
        self.model = model
        self.statements = statements

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
        # Make a copy of the model and add the observable for the modified
        # substrate
        new_model = deepcopy(self.model)
        sub_obs = Observable('target', sub_mp, _export=False)
        new_model.add_component(sub_obs)
        return True

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
            if embeds_into(rule_cp, cp):
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

