from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
import logging
import networkx
import itertools
import numpy as np
from copy import deepcopy
from collections import deque
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

    def add_statements(self, stmts):
        """Add to the list of statements to check against the model."""
        self.statements += stmts

    def get_im(self, force_update=True):
        """Get the influence map for the model, generating it if necessary.

        Parameters
        ----------
        force_update : bool
            Whether to generate the influence map when the function is called.
            If False, returns the previously generated influence map if
            available. Defaults to True.

        Returns
        -------
        pygraphviz AGraph object containing the influence map.
            The influence map can be rendered as a pdf using the dot layout
            program as follows::

                influence_map.draw('influence_map.pdf', prog='dot')
        """
        if self._im and not force_update:
            return self._im
        if not self.model:
            raise Exception("Cannot get influence map if there is no model.")
        else:
            logger.info("Generating influence map")
            self._im = kappa.influence_map(self.model)
            self._im.is_multigraph = lambda: False
            return self._im

    def check_model(self):
        """Check all the statements added to the ModelChecker.

        More efficient than check_statement when checking multiple statements
        because all relevant observables are added before building the
        influence map, preventing it from being repeatedly generated.

        Returns
        -------
        list of (Statement, bool)
            Each tuple contains the Statement checked against the model and
            a boolean value indicating whether the model can satisfies it.
        """
        results = []
        for stmt in self.statements:
            result = self.check_statement(stmt)
            results.append((stmt, result))
        return results

    def check_statement(self, stmt):
        """Check a single Statement against the model.

        Parameters
        ----------
        indra.statements.Statement
            The Statement to check.

        Returns
        -------
        boolean
            Whether the model satisfies the statement.
        """
        if isinstance(stmt, Modification):
            return self._check_modification(stmt)
        elif isinstance(stmt, RegulateActivity):
            return self._check_regulate_activity(stmt)
        else:
            return False

    def _check_regulate_activity(self, stmt):
        """Check a RegulateActivity statement."""
        logger.info('Checking stmt: %s' % stmt)
        # FIXME Currently this will match rules with the corresponding monomer
        # pattern from the Activation/Inhibition statement, which will nearly
        # always have no state conditions on it. In future, this statement foo
        # should also match rules in which 1) the agent is in its active form,
        # or 2) the agent is tagged as the enzyme in a rule of the appropriate
        # activity (e.g., a phosphorylation rule) FIXME
        subj_mp = pa.get_monomer_pattern(self.model, stmt.subj)

        target_polarity = 1 if stmt.is_activation else -1
        # This may fail, since there may be no rule in the model activating the
        # object, and the object may not have an "active" site of the
        # appropriate type
        obj_obs_name = pa.get_agent_rule_str(stmt.obj) + '_obs'
        try:
            obj_site_pattern = pa.get_site_pattern(stmt.obj)
            obj_site_pattern.update({stmt.obj_activity: 'active'})
            obj_monomer = self.model.monomers[stmt.obj.name]
            obj_mp = obj_monomer(**obj_site_pattern)
        except Exception as e:
            logger.info("Could not create obj monomer pattern: %s" % e)
            return False
        obj_obs = Observable(obj_obs_name, obj_mp, _export=False)
        return self._find_im_paths(subj_mp, obj_obs, target_polarity)

    def _check_modification(self, stmt):
        """Check a Modification statement."""
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
        mod_condition_name = modclass_to_modtype[stmt.__class__]
        if isinstance(stmt, RemoveModification):
            mod_condition_name = modtype_to_inverse[mod_condition_name]
        modified_sub = _add_modification_to_agent(stmt.sub, mod_condition_name,
                                                  stmt.residue, stmt.position)
        obs_name = pa.get_agent_rule_str(modified_sub) + '_obs'
        obj_mps = list(pa.grounded_monomer_patterns(self.model, modified_sub))
        if not obj_mps:
            logger.info('Failed to create observable; returning False')
            return False
        # Try to find paths between pairs of matching subj and object monomer
        # patterns
        for enz_mp, obj_mp in itertools.product(enz_mps, obj_mps):
            obj_obs = Observable(obs_name, obj_mp, _export=False)
            # Return True for the first valid path we find
            if self._find_im_paths(enz_mp, obj_obs, target_polarity):
                return True
        # If we got here, then there was no path for any observable
        return False

    def _find_im_paths(self, subj_mp, obj_obs, target_polarity):
        """Check for a source/target path in the influence map.

        Parameters
        ----------
        subj_mp : pysb.MonomerPattern
            MonomerPattern corresponding to the subject of the Statement
            being checked.
        obj_obs : pysb.Observable
            Observable corresponding to the object/target of the Statement
            being checked.
        target_polarity : 1 or -1
            Whether the influence in the Statement is positive (1) or negative
            (-1).

        Returns
        -------
        boolean
            Whether there is a path from a rule matching the subject
            MonomerPattern to the object Observable with the appropriate
            polarity.
        """
        # Reset the observables list
        self.model.observables = ComponentSet([])
        self.model.add_component(obj_obs)
        # Find rules in the model corresponding to the input
        logger.info('Finding paths between %s and %s with polarity %s' %
                    (subj_mp, obj_obs, target_polarity))
        if subj_mp is None:
            input_rule_set = None
        else:
            input_rules = _match_lhs(subj_mp, self.model.rules)
            logger.info('Found %s input rules matching %s' %
                        (len(input_rules), str(subj_mp)))
            # Filter to include only rules where the subj_mp is actually the
            # subject (i.e., don't pick up upstream rules where the subject
            # is itself a substrate/object)
            # FIXME: Note that this will eliminate rules where the subject
            # being checked is included on the left hand side as 
            # a bound condition rather than as an enzyme.
            subj_rules = pa.rules_with_annotation(self.model,
                                                  subj_mp.monomer.name,
                                                  'rule_has_subject')
            logger.info('%d rules with %s as subject' %
                        (len(subj_rules), subj_mp.monomer.name))
            input_rule_set = set([r.name for r in input_rules]).intersection(
                                 set([r.name for r in subj_rules]))
            logger.info('Final input rule set contains %d rules' %
                        len(input_rule_set))
            # If we have enzyme information but there are no input rules
            # matching the enzyme, then there is no path
            if not input_rule_set:
                return False
        # Generate the predecessors to our observable and count the paths
        # TODO: Make it optionally possible to return on the first path?
        num_paths = 0
        for path in _find_sources(self.get_im(), obj_obs.name, input_rule_set,
                                  target_polarity):
            num_paths += 1
        #for path in _find_sources_with_paths(self.get_im(), obj_obs.name,
        #                                     input_rule_set, target_polarity):
        #    num_paths += 1
        if num_paths > 0:
            return True
        else:
            return False

def _find_sources_with_paths(im, target, sources, polarity):
    """Get the subset of source nodes with paths to the target.

    Given a target, a list of sources, and a path polarity, perform a
    breadth-first search upstream from the target to find paths to any of the
    upstream sources.

    Parameters
    ----------
    im : pygraphviz.AGraph
        Graph containing the influence map.
    target : string
        The node (rule name) in the influence map to start looking upstream for
        marching sources.
    sources : list of strings
        The nodes (rules) corresponding to the subject or upstream influence
        being checked.
    polarity : int
        Required polarity of the path between source and target.

    Returns
    -------
    generator of path
        Yields paths as lists of nodes (rule names).  If there are no paths
        to any of the given source nodes, the generator is empty.
    """
    # First, create a list of visited nodes
    # Adapted from
    # http://stackoverflow.com/questions/8922060/
    #                       how-to-trace-the-path-in-a-breadth-first-search
    queue = deque([[target]])
    while queue:
        # Get the first path in the queue
        path = queue.popleft()
        node = path[-1]
        # If there's only one node in the path, it's the observable we're
        # starting from, so the path is positive
        if len(path) == 1:
            sign = 1
        # Because the path runs from target back to source, we have to reverse
        # the path to calculate the overall polarity
        else:
            sign = _path_polarity(im, reversed(path))
        if (sources is None or node in sources) and sign == polarity:
            logger.info('Found path: %s' % path)
            yield path
        for predecessor, sign in _get_signed_predecessors(im, node, 1):
            new_path = list(path)
            new_path.append(predecessor)
            queue.append(new_path)
    return

def _find_sources(im, target, sources, polarity):
    """Get the subset of source nodes with paths to the target.

    Given a target, a list of sources, and a path polarity, perform a
    breadth-first search upstream from the target to determine whether any of
    the queried sources have paths to the target with the appropriate polarity.
    For efficiency, does not return the full path, but identifies the upstream
    sources and the length of the path.

    Parameters
    ----------
    im : pygraphviz.AGraph
        Graph containing the influence map.
    target : string
        The node (rule name) in the influence map to start looking upstream for
        marching sources.
    sources : list of strings
        The nodes (rules) corresponding to the subject or upstream influence
        being checked.
    polarity : int
        Required polarity of the path between source and target.

    Returns
    -------
    generator of (source, polarity, path_length)
        Yields tuples of source node (string), polarity (int) and path length
        (int). If there are no paths to any of the given source nodes, the
        generator isignempty.
    """
    # First, create a list of visited nodes
    # Adapted from
    # networkx.algorithms.traversal.breadth_first_search.bfs_edges
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
            # Is this child one of the source nodes we're looking for? If so,
            # yield it along with path length
            if (sources is None or child in sources) and sign == polarity:
                logger.info("Found path to %s from %s with desired sign %s "
                            "with length %d" %
                            (target, child, polarity, path_length))
                yield (child, sign, path_length)
            # Check this child against the visited list. If we haven't visited
            # it already (accounting for the path to the node), then add it
            # to the queue.
            if (child, sign) not in visited:
                visited.add((child, sign))
                queue.append((child, _get_signed_predecessors(im, child, sign),
                              path_length + 1))
        # Once we've finished iterating over the children of the current node,
        # pop the node off and go to the next one in the queue
        except StopIteration:
            queue.popleft()
            path_length += 1
    # There was no path; this will produce an empty generator
    return


def _get_signed_predecessors(im, node, polarity):
    """Get upstream nodes in the influence map

    Return the upstream nodes along with the overall polarity of the path
    to that node by account for the polarity of the path to the given node
    and the polarity of the edge between the given node and its immediate
    predecessors.

    Parameters
    ----------
    im : pygraphviz.AGraph
        Graph containing the influence map.
    node : string
        The node (rule name) in the influence map to get predecessors (upstream
        nodes) for.
    polarity : int
        Polarity of the overall path to the given node.


    Returns
    -------
    generator of tuples, (node, polarity)
        Each tuple returned contains two elements, a node (string) and the
        polarity of the overall path (int) to that node.
    """
    signed_pred_list = []
    predecessors = im.predecessors_iter
    for pred in predecessors(node):
        pred_edge = im.get_edge(pred, node)
        yield (pred, _get_edge_sign(pred_edge) * polarity)


def _get_edge_sign(edge):
    """Get the polarity of the influence by examining the edge color."""
    if edge.attr.get('color') is None:
        raise Exception('No color attribute for edge.')
    elif edge.attr['color'] == 'green':
        return 1
    elif edge.attr['color'] == 'red':
        return -1
    else:
        raise Exception('Unexpected edge color: %s' % edge.attr['color'])


def _add_modification_to_agent(agent, mod_type, residue, position):
    """Add a modification condition to an Agent."""
    new_mod = ModCondition(mod_type, residue, position)
    # Check if this modification already exists
    for old_mod in agent.mods:
        if old_mod.equals(new_mod):
            return agent
    new_agent = deepcopy(agent)
    new_agent.mods.append(new_mod)
    return new_agent


def _match_lhs(cp, rules):
    """Get rules with a left-hand side matching the given ComplexPattern."""
    rule_matches = []
    for rule in rules:
        reactant_pattern = rule.rule_expression.reactant_pattern
        for rule_cp in reactant_pattern.complex_patterns:
            if _cp_embeds_into(rule_cp, cp):
                rule_matches.append(rule)
                break
    return rule_matches


def _cp_embeds_into(cp1, cp2):
    """Check that any state in ComplexPattern2 is matched in ComplexPattern1.
    """
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
            if _mp_embeds_into(mp1, mp2):
                return True
    return False


def _mp_embeds_into(mp1, mp2):
    """Check that conditions in MonomerPattern2 are met in MonomerPattern1."""
    sc_matches = []
    if mp1.monomer.name != mp2.monomer.name:
        return False
    # Check that all conditions in mp2 are met in mp1
    for site_name, site_state in mp2.site_conditions.items():
        if site_name not in mp1.site_conditions or \
           site_state != mp1.site_conditions[site_name]:
            return False
    return True


"""
# NOTE: This code is currently "deprecated" because it has been replaced by the
# use of Observables for the Statement objects.

def match_rhs(cp, rules):
    rule_matches = []
    for rule in rules:
        product_pattern = rule.rule_expression.product_pattern
        for rule_cp in product_pattern.complex_patterns:
            if _cp_embeds_into(rule_cp, cp):
                rule_matches.append(rule)
                break
    return rule_matches

def find_production_rules(cp, rules):
    # Find rules where the CP matches the left hand side
    lhs_rule_set = set(_match_lhs(cp, rules))
    # Now find rules where the CP matches the right hand side
    rhs_rule_set = set(match_rhs(cp, rules))
    # Production rules are rules where there is a match on the right hand
    # side but not on the left hand side
    prod_rules = list(rhs_rule_set.difference(lhs_rule_set))
    return prod_rules

def find_consumption_rules(cp, rules):
    # Find rules where the CP matches the left hand side
    lhs_rule_set = set(_match_lhs(cp, rules))
    # Now find rules where the CP matches the right hand side
    rhs_rule_set = set(match_rhs(cp, rules))
    # Consumption rules are rules where there is a match on the left hand
    # side but not on the right hand side
    cons_rules = list(lhs_rule_set.difference(rhs_rule_set))
    return cons_rules
"""
def _path_polarity(im, path):
    # This doesn't address the effect of the rules themselves on the
    # observables of interest--just the effects of the rules on each other
    edge_polarities = []
    path_list = list(path)
    edges = zip(path_list[0:-1], path_list[1:])
    for from_rule, to_rule in edges:
        edge = im.get_edge(from_rule, to_rule)
        edge_polarities.append(_get_edge_sign(edge))
    # Compute and return the overall path polarity
    path_polarity = np.prod(edge_polarities)
    assert path_polarity == 1 or path_polarity == -1
    #return True if path_polarity == 1 else False
    return path_polarity

