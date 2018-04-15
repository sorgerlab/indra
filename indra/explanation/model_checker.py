from __future__ import print_function, unicode_literals, absolute_import
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import logging
import numbers
import textwrap
import networkx as nx
import itertools
import numpy as np
import scipy.stats
from copy import deepcopy
from collections import deque, defaultdict, namedtuple
import kappy
from pysb import WILD, export, Observable, ComponentSet
from pysb.core import as_complex_pattern, ComponentDuplicateNameError
from indra.statements import *
from indra.assemblers import pysb_assembler as pa
from indra.tools.expand_families import _agent_from_uri
from indra.explanation import paths_graph as pg
from collections import Counter
from indra.util.kappa_util import im_json_to_graph

logger = logging.getLogger('model_checker')


class PathMetric(object):
    """Describes results of simple path search (path existence)."""
    def __init__(self, source_node, target_node, polarity, length):
        self.source_node = source_node
        self.target_node = target_node
        self.polarity = polarity
        self.length = length

    def __repr__(self):
        return str(self)

    @python_2_unicode_compatible
    def __str__(self):
        return ('source_node: %s, target_node: %s, polarity: %s, length: %d' %
                (self.source_node, self.target_node, self.polarity,
                 self.length))

class PathResult(object):
    """Describes results of running the ModelChecker on a single Statement.

    Parameters
    ----------
    path_found : bool
    result_code : string
        STATEMENT_TYPE_NOT_HANDLED
        SUBJECT_MONOMERS_NOT_FOUND
        OBSERVABLES_NOT_FOUND
        NO_PATHS_FOUND
        MAX_PATH_LENGTH_EXCEEDED
        PATHS_FOUND
        INPUT_RULES_NOT_FOUND
        MAX_PATHS_ZERO

    Attributes
    ----------
    path_found : boolean
    result_code : string
    path_metrics : list of PathMetric
    paths : list of paths
    max_paths :
    max_path_length :
    """
    def __init__(self, path_found, result_code, max_paths, max_path_length):
        self.path_found = path_found
        self.result_code = result_code
        self.max_paths = max_paths
        self.max_path_length = max_path_length
        self.path_metrics = []
        self.paths = []

    def add_path(self, path):
        self.paths.append(path)

    def add_metric(self, path_metric):
        self.path_metrics.append(path_metric)

    @python_2_unicode_compatible
    def __str__(self):
        summary = textwrap.dedent("""
            PathResult:
                path_found: {path_found}
                result_code: {result_code}
                path_metrics: {path_metrics}
                paths: {paths}
                max_paths: {max_paths}
                max_path_length: {max_path_length}""")
        ws = '\n        '
        # String representation of path metrics
        if not self.path_metrics:
            pm_str = str(self.path_metrics)
        else:
            pm_str = ws + ws.join(['%d: %s' % (pm_ix, pm) for pm_ix, pm in
                                            enumerate(self.path_metrics)])
        def format_path(path, num_spaces=11):
            path_ws = '\n' + (' ' * num_spaces)
            return path_ws.join([str(p) for p in path])

        # String representation of paths
        if not self.paths:
            path_str = str(self.paths)
        else:
            path_str = ws + ws.join(['%d: %s' % (p_ix, format_path(p))
                                     for p_ix, p in enumerate(self.paths)])

        return summary.format(path_found=self.path_found,
                       result_code=self.result_code,
                       max_paths=self.max_paths,
                       max_path_length=self.max_path_length,
                       path_metrics=pm_str, paths=path_str)

    def __repr__(self):
        return str(self)

class ModelChecker(object):
    """Check a PySB model against a set of INDRA statements.

    Parameters
    ----------
    model : pysb.Model
        A PySB model to check.
    statements : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to check the model against.
    agent_obs: Optional[list[indra.statements.Agent]]
        A list of INDRA Agents in a given state to be observed.
    do_sampling : bool
        Whether to use breadth-first search or weighted sampling to
        generate paths. Default is False (breadth-first search).
    seed : int
        Random seed for sampling (optional, default is None).
    """

    def __init__(self, model, statements=None, agent_obs=None,
                 do_sampling=False, seed=None):
        self.model = model
        if statements:
            self.statements = statements
        else:
            self.statements = []
        if agent_obs:
            self.agent_obs = agent_obs
        else:
            self.agent_obs = []
        if seed is not None:
            np.random.seed(seed)
        # Whether to do sampling
        self.do_sampling = do_sampling
        # Influence map
        self._im = None
        # Map from statements to associated observables
        self.stmt_to_obs = {}
        # Map from agents to associated observables
        self.agent_to_obs = {}
        # Map between rules and downstream observables
        self.rule_obs_dict = {}

    def add_statements(self, stmts):
        """Add to the list of statements to check against the model.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            The list of Statements to be added for checking.
        """
        self.statements += stmts

    def generate_im(self, model):
        """Return a graph representing the influence map generated by Kappa

        Parameters
        ----------
        model : pysb.Model
            The PySB model whose influence map is to be generated

        Returns
        -------
        graph : networkx.MultiDiGraph
            A MultiDiGraph representing the influence map
        """
        kappa = kappy.KappaStd()
        model_str = export.export(model, 'kappa')
        kappa.add_model_string(model_str)
        kappa.project_parse()
        imap = kappa.analyses_influence_map()
        graph = im_json_to_graph(imap)
        return graph

    def get_im(self, force_update=False):
        """Get the influence map for the model, generating it if necessary.

        Parameters
        ----------
        force_update : bool
            Whether to generate the influence map when the function is called.
            If False, returns the previously generated influence map if
            available. Defaults to True.

        Returns
        -------
        networkx MultiDiGraph object containing the influence map.
            The influence map can be rendered as a pdf using the dot layout
            program as follows::

                im_agraph = nx.nx_agraph.to_agraph(influence_map)
                im_agraph.draw('influence_map.pdf', prog='dot')
        """
        if self._im and not force_update:
            return self._im
        if not self.model:
            raise Exception("Cannot get influence map if there is no model.")

        def add_obs_for_agent(agent):
            obj_mps = list(pa.grounded_monomer_patterns(self.model, agent))
            if not obj_mps:
                logger.debug('No monomer patterns found in model for agent %s, '
                            'skipping' % agent)
                return
            obs_list = []
            for obj_mp in obj_mps:
                obs_name = _monomer_pattern_label(obj_mp) + '_obs'
                # Add the observable
                obj_obs = Observable(obs_name, obj_mp, _export=False)
                obs_list.append(obs_name)
                try:
                    self.model.add_component(obj_obs)
                except ComponentDuplicateNameError as e:
                    pass
            return obs_list

        # Create observables for all statements to check, and add to model
        # Remove any existing observables in the model
        self.model.observables = ComponentSet([])
        for stmt in self.statements:
            # Generate observables for Modification statements
            if isinstance(stmt, Modification):
                mod_condition_name = modclass_to_modtype[stmt.__class__]
                if isinstance(stmt, RemoveModification):
                    mod_condition_name = modtype_to_inverse[mod_condition_name]
                # Add modification to substrate agent
                modified_sub = _add_modification_to_agent(stmt.sub,
                                    mod_condition_name, stmt.residue,
                                    stmt.position)
                obs_list = add_obs_for_agent(modified_sub)
                # Associate this statement with this observable
                self.stmt_to_obs[stmt] = obs_list
            # Generate observables for Activation/Inhibition statements
            elif isinstance(stmt, RegulateActivity):
                regulated_obj, polarity = \
                        _add_activity_to_agent(stmt.obj, stmt.obj_activity,
                                               stmt.is_activation)
                obs_list = add_obs_for_agent(regulated_obj)
                # Associate this statement with this observable
                self.stmt_to_obs[stmt] = obs_list
            elif isinstance(stmt, RegulateAmount):
                obs_list = add_obs_for_agent(stmt.obj)
                self.stmt_to_obs[stmt] = obs_list
        # Add observables for each agent
        for ag in self.agent_obs:
            obs_list = add_obs_for_agent(ag)
            self.agent_to_obs[ag] = obs_list

        logger.info("Generating influence map")
        self._im = self.generate_im(self.model)
        #self._im.is_multigraph = lambda: False
        # Now, for every rule in the model, check if there are any observables
        # downstream; alternatively, for every observable in the model, get a
        # list of rules.
        # We'll need the dictionary to check if nodes are observables
        node_attributes = nx.get_node_attributes(self._im, 'node_type')
        for rule in self.model.rules:
            obs_list = []
            # Get successors of the rule node
            for neighb in self._im.neighbors(rule.name):
                # Check if the node is an observable
                if node_attributes[neighb] != 'variable':
                    continue
                # Get the edge and check the polarity
                edge_sign = _get_edge_sign(self._im, (rule.name, neighb))
                obs_list.append((neighb, edge_sign))
            self.rule_obs_dict[rule.name] = obs_list
        return self._im

    def check_model(self, max_paths=1, max_path_length=5):
        """Check all the statements added to the ModelChecker.

        Parameters
        ----------
        max_paths : Optional[int]
            The maximum number of specific paths to return for each Statement
            to be explained. Default: 1
        max_path_length : Optional[int]
            The maximum length of specific paths to return. Default: 5

        Returns
        -------
        list of (Statement, PathResult)
            Each tuple contains the Statement checked against the model and
            a PathResult object describing the results of model checking.
        """
        results = []
        for stmt in self.statements:
            result = self.check_statement(stmt, max_paths, max_path_length)
            results.append((stmt, result))
        return results

    def check_statement(self, stmt, max_paths=1, max_path_length=5):
        """Check a single Statement against the model.

        Parameters
        ----------
        stmt : indra.statements.Statement
            The Statement to check.
        max_paths : Optional[int]
            The maximum number of specific paths to return for each Statement
            to be explained. Default: 1
        max_path_length : Optional[int]
            The maximum length of specific paths to return. Default: 5

        Returns
        -------
        boolean
            True if the model satisfies the Statement.
        """
        # Make sure the influence map is initialized
        self.get_im()
        if isinstance(stmt, Modification):
            return self._check_modification(stmt, max_paths, max_path_length)
        elif isinstance(stmt, RegulateActivity):
            return self._check_regulate_activity(stmt, max_paths,
                                                 max_path_length)
        elif isinstance(stmt, RegulateAmount):
            return self._check_regulate_amount(stmt, max_paths,
                                               max_path_length)
        else:
            return PathResult(False, 'STATEMENT_TYPE_NOT_HANDLED',
                              max_paths, max_path_length)

    def _check_regulate_activity(self, stmt, max_paths, max_path_length):
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
        obs_names = self.stmt_to_obs[stmt]
        for obs_name in obs_names:
            return self._find_im_paths(subj_mp, obs_name, target_polarity,
                                       max_paths, max_path_length)

    def _check_regulate_amount(self, stmt, max_paths, max_path_length):
        """Check a RegulateAmount statement."""
        logger.info('Checking stmt: %s' % stmt)
        subj_mp = pa.get_monomer_pattern(self.model, stmt.subj)
        target_polarity = 1 if isinstance(stmt, IncreaseAmount) else -1
        obs_names = self.stmt_to_obs[stmt]
        for obs_name in obs_names:
            return self._find_im_paths(subj_mp, obs_name, target_polarity,
                                       max_paths, max_path_length)

    def _check_modification(self, stmt, max_paths, max_path_length):
        """Check a Modification statement."""
        # Identify the observable we're looking for in the model, which
        # may not exist!
        # The observable is the modified form of the substrate
        logger.info('Checking stmt: %s' % stmt)
        # Look for an agent with the appropriate grounding in the model
        if stmt.enz is not None:
            enz_mps = list(pa.grounded_monomer_patterns(self.model, stmt.enz))
            if not enz_mps:
                logger.debug('No monomers found corresponding to agent %s' %
                             stmt.enz)
                return PathResult(False, 'SUBJECT_MONOMERS_NOT_FOUND',
                                  max_paths, max_path_length)
        else:
            enz_mps = [None]
        # Get target polarity
        target_polarity = -1 if isinstance(stmt, RemoveModification) else 1
        obs_names = self.stmt_to_obs[stmt]
        if not obs_names:
            logger.debug("No observables for stmt %s, returning False" % stmt)
            return PathResult(False, 'OBSERVABLES_NOT_FOUND',
                              max_paths, max_path_length)

        for enz_mp, obs_name in itertools.product(enz_mps, obs_names):
            # FIXME Returns on the path found for the first enz_mp/obs combo
            result = self._find_im_paths(enz_mp, obs_name, target_polarity,
                                         max_paths, max_path_length)
            # If result for this observable is not False, then we return it;
            # otherwise, that means there was no path for this observable, so
            # we have to try the next one
            if result.path_found:
                return result
        # If we got here, then there was no path for any observable
        return PathResult(False, 'NO_PATHS_FOUND',
                          max_paths, max_path_length)

    def _get_input_rules(self, subj_mp):
        if subj_mp is None:
            raise ValueError("Cannot take None as an argument for subj_mp.")
        input_rules = _match_lhs(subj_mp, self.model.rules)
        logger.debug('Found %s input rules matching %s' %
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
        logger.debug('%d rules with %s as subject' %
                     (len(subj_rules), subj_mp.monomer.name))
        input_rule_set = set([r.name for r in input_rules]).intersection(
                             set([r.name for r in subj_rules]))
        logger.debug('Final input rule set contains %d rules' %
                     len(input_rule_set))
        return input_rule_set

    def _sample_paths(self, input_rule_set, obs_name, target_polarity,
                      max_paths=1, max_path_length=5):
        if max_paths == 0:
            raise ValueError("max_paths cannot be 0 for path sampling.")
        # Convert path polarity representation from 0/1 to 1/-1
        def convert_polarities(path_list):
            return [tuple((n[0], 0 if n[1] > 0 else 1)
                          for n in path)
                          for path in path_list]

        pg_polarity = 0 if target_polarity > 0 else 1
        nx_graph = _im_to_signed_digraph(self.get_im())
        # Add edges from dummy node to input rules
        source_node = 'SOURCE_NODE'
        for rule in input_rule_set:
            nx_graph.add_edge(source_node, rule, attr_dict={'sign': 0})
        # -------------------------------------------------
        # Create combined paths_graph
        f_level, b_level = pg.get_reachable_sets(nx_graph, source_node,
                                                 obs_name, max_path_length,
                                                 signed=True)
        pg_list = []
        for path_length in range(1, max_path_length+1):
            cfpg = pg.CFPG.from_graph(
                    nx_graph, source_node, obs_name, path_length, f_level,
                    b_level, signed=True, target_polarity=pg_polarity)
            pg_list.append(cfpg)
        combined_pg = pg.CombinedCFPG(pg_list)
        # Make sure the combined paths graph is not empty
        if not combined_pg.graph:
            pr = PathResult(False, 'NO_PATHS_FOUND', max_paths, max_path_length)
            pr.path_metrics = None
            pr.paths = []
            return pr

        # Get a dict of rule objects
        rule_obj_dict = {}
        for ann in self.model.annotations:
            if ann.predicate == 'rule_has_object':
                rule_obj_dict[ann.subject] = ann.object

        # Get monomer initial conditions
        ic_dict = {}
        for mon in self.model.monomers:
            # FIXME: A hack that depends on the _0 convention
            ic_name = '%s_0' % mon.name
            # TODO: Wrap this in try/except?
            ic_param = self.model.parameters[ic_name]
            ic_value = ic_param.value
            ic_dict[mon.name] = ic_value

        # Set weights in PG based on model initial conditions
        for cur_node in combined_pg.graph.nodes():
            edge_weights = {}
            rule_obj_list = []
            edge_weights_by_gene = {}
            for u, v in combined_pg.graph.out_edges(cur_node):
                v_rule = v[1][0]
                # Get the object of the rule (a monomer name)
                rule_obj = rule_obj_dict.get(v_rule)
                if rule_obj:
                    # Add to list so we can count instances by gene
                    rule_obj_list.append(rule_obj)
                    # Get the abundance of rule object from the initial
                    # conditions
                    # TODO: Wrap in try/except?
                    ic_value = ic_dict[rule_obj]
                else:
                    ic_value = 1.0
                edge_weights[(u, v)] = ic_value
                edge_weights_by_gene[rule_obj] = ic_value
            # Get frequency of different rule objects
            rule_obj_ctr = Counter(rule_obj_list)
            # Normalize results by weight sum and gene frequency at this level
            edge_weight_sum = sum(edge_weights_by_gene.values())
            edge_weights_norm = {}
            for e, v in edge_weights.items():
                v_rule = e[1][1][0]
                rule_obj = rule_obj_dict.get(v_rule)
                if rule_obj:
                    rule_obj_count = rule_obj_ctr[rule_obj]
                else:
                    rule_obj_count = 1
                edge_weights_norm[e] = ((v / float(edge_weight_sum)) /
                                        float(rule_obj_count))
            # Add edge weights to paths graph
            nx.set_edge_attributes(combined_pg.graph, 'weight',
                                   edge_weights_norm)

        # Sample from the combined CFPG
        paths = combined_pg.sample_paths(max_paths)
        # -------------------------------------------------
        if paths:
            pr = PathResult(True, 'PATHS_FOUND', max_paths, max_path_length)
            pr.path_metrics = None
            # Convert path polarity representation from 0/1 to 1/-1
            pr.paths = convert_polarities(paths)
            # Strip off the SOURCE_NODE prefix
            pr.paths = [p[1:] for p in pr.paths]
        else:
            assert False
            pr = PathResult(False, 'NO_PATHS_FOUND', max_paths, max_path_length)
            pr.path_metrics = None
            pr.paths = []
        return pr

    def _find_im_paths(self, subj_mp, obs_name, target_polarity,
                       max_paths=1, max_path_length=5):
        """Check for a source/target path in the influence map.

        Parameters
        ----------
        subj_mp : pysb.MonomerPattern
            MonomerPattern corresponding to the subject of the Statement
            being checked.
        obs_name : string
            Name of the PySB model Observable corresponding to the
            object/target of the Statement being checked.
        target_polarity : 1 or -1
            Whether the influence in the Statement is positive (1) or negative
            (-1).

        Returns
        -------
        boolean or list of str
            Whether there is a path from a rule matching the subject
            MonomerPattern to the object Observable with the appropriate
            polarity.
        """
        logger.info(('Running path finding with max_paths=%d,'
                     ' max_path_length=%d') % (max_paths, max_path_length))
        # Find rules in the model corresponding to the input
        if subj_mp is None:
            input_rule_set = None
        else:
            input_rule_set = self._get_input_rules(subj_mp)
            if not input_rule_set:
                return PathResult(False, 'INPUT_RULES_NOT_FOUND',
                                  max_paths, max_path_length)
        logger.info('Finding paths between %s and %s with polarity %s' %
                    (subj_mp, obs_name, target_polarity))

        # -- Route to the path sampling function --
        if self.do_sampling:
            return self._sample_paths(input_rule_set, obs_name, target_polarity,
                               max_paths, max_path_length)

        # -- Do Breadth-First Enumeration --
        # Generate the predecessors to our observable and count the paths
        path_lengths = []
        path_metrics = []
        for source, polarity, path_length in \
                    _find_sources(self.get_im(), obs_name, input_rule_set,
                                  target_polarity):

            pm = PathMetric(source, obs_name, polarity, path_length)
            path_metrics.append(pm)
            path_lengths.append(path_length)
        logger.info('Finding paths between %s and %s with polarity %s' %
                    (subj_mp, obs_name, target_polarity))
        # Now, look for paths
        paths = []
        if path_metrics and max_paths == 0:
            pr = PathResult(True, 'MAX_PATHS_ZERO',
                            max_paths, max_path_length)
            pr.path_metrics = path_metrics
            return pr
        elif path_metrics:
            if min(path_lengths) <= max_path_length:
                pr = PathResult(True, 'PATHS_FOUND', max_paths, max_path_length)
                pr.path_metrics = path_metrics
                # Get the first path
                path_iter = enumerate(_find_sources_with_paths(
                                           self.get_im(), obs_name,
                                           input_rule_set, target_polarity))
                for path_ix, path in path_iter:
                    flipped = _flip(self.get_im(), path)
                    pr.add_path(flipped)
                    if len(pr.paths) >= max_paths:
                        break
                return pr
            # There are no paths shorter than the max path length, so we
            # don't bother trying to get them
            else:
                pr = PathResult(True, 'MAX_PATH_LENGTH_EXCEEDED',
                                max_paths, max_path_length)
                pr.path_metrics = path_metrics
                return pr
        else:
            return PathResult(False, 'NO_PATHS_FOUND',
                              max_paths, max_path_length)

    def score_paths(self, paths, agents_values, loss_of_function=False,
                    sigma=0.15, include_final_node=False):
        """Return scores associated with a given set of paths.

        Parameters
        ----------
        paths : list[list[tuple[str, int]]]
            A list of paths obtained from path finding. Each path is a list
            of tuples (which are edges in the path), with the first element
            of the tuple the name of a rule, and the second element its
            polarity in the path.
        agents_values : dict[indra.statements.Agent, float]
            A dictionary of INDRA Agents and their corresponding measured
            value in a given experimental condition.
        loss_of_function : Optional[boolean]
            If True, flip the polarity of the path. For instance, if the effect
            of an inhibitory drug is explained, set this to True.
            Default: False
        sigma : Optional[float]
            The estimated standard deviation for the normally distributed
            measurement error in the observation model used to score paths
            with respect to data. Default: 0.15
        include_final_node : Optional[boolean]
            Determines whether the final node of the path is included in the
            score. Default: False
        """
        obs_model = lambda x: scipy.stats.norm(x, sigma)
        # Build up dict mapping observables to values
        obs_dict = {}
        for ag, val in agents_values.items():
            obs_list = self.agent_to_obs[ag]
            if obs_list is not None:
                for obs in obs_list:
                    obs_dict[obs] = val
        # For every path...
        path_scores = []
        for path in paths:
            logger.info('------')
            logger.info("Scoring path:")
            logger.info(path)
            # Look at every node in the path, excluding the final
            # observable...
            path_score = 0
            last_path_node_index = -1 if include_final_node else -2
            for node, sign in path[:last_path_node_index]:
                # ...and for each node check the sign to see if it matches the
                # data. So the first thing is to look at what's downstream
                # of the rule
                # affected_obs is a list of observable names alogn
                for affected_obs, rule_obs_sign in self.rule_obs_dict[node]:
                    flip_polarity = -1 if loss_of_function else 1
                    pred_sign = sign * rule_obs_sign * flip_polarity
                    # Check to see if this observable is in the data
                    logger.info('%s %s: effect %s %s' %
                                (node, sign, affected_obs, pred_sign))
                    measured_val = obs_dict.get(affected_obs)
                    if measured_val:
                        # For negative predictions use CDF (prob that given
                        # measured value, true value lies below 0)
                        if pred_sign <= 0:
                            prob_correct = obs_model(measured_val).logcdf(0)
                        # For positive predictions, use log survival function
                        # (SF = 1 - CDF, i.e., prob that true value is
                        # above 0)
                        else:
                            prob_correct = obs_model(measured_val).logsf(0)
                        logger.info('Actual: %s, Log Probability: %s' %
                                    (measured_val, prob_correct))
                        path_score += prob_correct
                if not self.rule_obs_dict[node]:
                    logger.info('%s %s' % (node, sign))
                    prob_correct = obs_model(0).logcdf(0)
                    logger.info('Unmeasured node, Log Probability: %s' %
                                (prob_correct))
                    path_score += prob_correct
            # Normalized path
            #path_score = path_score / len(path)
            logger.info("Path score: %s" % path_score)
            path_scores.append(path_score)
        path_tuples = list(zip(paths, path_scores))
        # Sort first by path length
        sorted_by_length = sorted(path_tuples, key=lambda x: len(x[0]))
        # Sort by probability; sort in reverse order to large values
        # (higher probabilities) are ranked higher
        scored_paths = sorted(sorted_by_length, key=lambda x: x[1],
                              reverse=True)
        return scored_paths

    def prune_influence_map(self):
        """Remove edges between rules causing problematic non-transitivity.

        First, all self-loops are removed. After this initial step, edges are
        removed between rules when they share *all* child nodes except for each
        other; that is, they have a mutual relationship with each other and
        share all of the same children.

        Note that edges must be removed in batch at the end to prevent edge
        removal from affecting the lists of rule children during the comparison
        process.
        """
        im = self.get_im()

        # First, remove all self-loops
        logger.info('Removing self loops')
        for e in im.edges():
            if e[0] == e[1]:
                logger.info('Removing self loop: %s', e)
                im.remove_edge(e[0], e[1])

        # Remove parameter nodes from influence map
        remove_im_params(self.model, im)

        # Now compare nodes pairwise and look for overlap between child nodes
        logger.info('Get successorts of each node')
        successors = im.successors_iter
        succ_dict = {}
        for node in im.nodes():
            succ_dict[node] = set(successors(node))
        # Sort and then group nodes by number of successors
        logger.info('Compare combinations of successors')
        group_key_fun = lambda x: len(succ_dict[x])
        nodes_sorted = sorted(im.nodes(), key=group_key_fun)
        groups = itertools.groupby(nodes_sorted, key=group_key_fun)
        # Now iterate over each group and then construct combinations
        # within the group to check for shared sucessors
        edges_to_remove = []
        for gix, group in groups:
            combos = itertools.combinations(group, 2)
            for ix, (p1, p2) in enumerate(combos):
                # Children are identical except for mutual relationship
                if succ_dict[p1].difference(succ_dict[p2]) == set([p2]) and \
                   succ_dict[p2].difference(succ_dict[p1]) == set([p1]):
                    for u, v in ((p1, p2), (p2, p1)):
                        edges_to_remove.append((u, v))
                        logger.debug('Will remove edge (%s, %s)', u, v)
        logger.info('Removing %d edges from influence map' %
                    len(edges_to_remove))
        # Now remove all the edges to be removed with a single call
        im.remove_edges_from(edges_to_remove)


def _find_sources_sample(im, target, sources, polarity, rule_obs_dict,
                         agent_to_obs, agents_values):
    # Build up dict mapping observables to values
    obs_dict = {}
    for ag, val in agents_values.items():
        obs_list = agent_to_obs[ag]
        for obs in obs_list:
            obs_dict[obs] = val

    sigma = 0.2
    def obs_model(x):
        return scipy.stats.norm(x, sigma)

    def _sample_pred(im, target, rule_obs_dict, obs_model):
        preds = list(_get_signed_predecessors(im, target, 1))
        if not preds:
            return None
        pred_scores = []
        for pred, sign in preds:
            pred_score = 0
            for affected_obs, rule_obs_sign in rule_obs_dict[pred]:
                pred_sign = sign * rule_obs_sign
                # Check to see if this observable is in the data
                logger.info('%s %s: effect %s %s' %
                            (pred, sign, affected_obs, pred_sign))
                measured_val = obs_dict.get(affected_obs)
                if measured_val:
                    logger.info('Actual: %s' % measured_val)
                    # The tail probability of the real value being above 1
                    tail_prob = obs_model(measured_val).cdf(1)
                    pred_score += (tail_prob if pred_sign == 1 else
                                   1-tail_prob)
            pred_scores.append(pred_score)
        # Normalize scores
        pred_scores = np.array(pred_scores) / np.sum(pred_scores)
        pred_idx = np.random.choice(range(len(preds)), p=pred_scores)
        pred = preds[pred_idx]
        return pred

    preds = []
    for i in range(100):
        pred = _sample_pred(im, target, rule_obs_dict, obs_model)
        preds.append(pred[0])

def _find_sources_with_paths(im, target, sources, polarity):
    """Get the subset of source nodes with paths to the target.

    Given a target, a list of sources, and a path polarity, perform a
    breadth-first search upstream from the target to find paths to any of the
    upstream sources.

    Parameters
    ----------
    im : networkx.MultiDiGraph
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
    # FIXME: the sign information for the target should be associated with
    # the observable itself
    queue = deque([[(target, 1)]])
    while queue:
        # Get the first path in the queue
        path = queue.popleft()
        node, node_sign = path[-1]
        # If there's only one node in the path, it's the observable we're
        # starting from, so the path is positive
        # if len(path) == 1:
        #    sign = 1
        # Because the path runs from target back to source, we have to reverse
        # the path to calculate the overall polarity
        #else:
        #    sign = _path_polarity(im, reversed(path))
        # Don't allow trivial paths consisting only of the target observable
        if (sources is None or node in sources) and node_sign == polarity \
           and len(path) > 1:
            logger.debug('Found path: %s' % str(_flip(im, path)))
            yield tuple(path)
        for predecessor, sign in _get_signed_predecessors(im, node, node_sign):
            # Only add predecessors to the path if it's not already in the
            # path--prevents loops
            if (predecessor, sign) in path:
                continue
            # Otherwise, the new path is a copy of the old one plus the new
            # predecessor
            new_path = list(path)
            new_path.append((predecessor, sign))
            queue.append(new_path)
    return


def remove_im_params(model, im):
    """Remove parameter nodes from the influence map.

    Parameters
    ----------
    model : pysb.core.Model
        PySB model.
    im : networkx.MultiDiGraph
        Influence map.

    Returns
    -------
    networkx.MultiDiGraph
        Influence map with the parameter nodes removed.
    """
    for param in model.parameters:
        # If the node doesn't exist e.g., it may have already been removed),
        # skip over the parameter without error
        try:
            im.remove_node(param.name)
        except:
            pass

def _find_sources(im, target, sources, polarity):
    """Get the subset of source nodes with paths to the target.

    Given a target, a list of sources, and a path polarity, perform a
    breadth-first search upstream from the target to determine whether any of
    the queried sources have paths to the target with the appropriate polarity.
    For efficiency, does not return the full path, but identifies the upstream
    sources and the length of the path.

    Parameters
    ----------
    im : networkx.MultiDiGraph
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
    queue = deque([(target_tuple, _get_signed_predecessors(im, target, 1), 0)])
    while queue:
        parent, children, path_length = queue[0]
        try:
            # Get the next child in the list
            (child, sign) = next(children)
            # Is this child one of the source nodes we're looking for? If so,
            # yield it along with path length.
            if (sources is None or child in sources) and sign == polarity:
                logger.debug("Found path to %s from %s with desired sign %s "
                             "with length %d" %
                             (target, child, polarity, path_length+1))
                yield (child, sign, path_length+1)
            # Check this child against the visited list. If we haven't visited
            # it already (accounting for the path to the node), then add it
            # to the queue.
            if (child, sign) not in visited:
                visited.add((child, sign))
                queue.append(((child, sign),
                              _get_signed_predecessors(im, child, sign),
                              path_length + 1))
        # Once we've finished iterating over the children of the current node,
        # pop the node off and go to the next one in the queue
        except StopIteration:
            queue.popleft()
    # There was no path; this will produce an empty generator
    return


def _get_signed_predecessors(im, node, polarity):
    """Get upstream nodes in the influence map.

    Return the upstream nodes along with the overall polarity of the path
    to that node by account for the polarity of the path to the given node
    and the polarity of the edge between the given node and its immediate
    predecessors.

    Parameters
    ----------
    im : networkx.MultiDiGraph
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
        pred_edge = (pred, node)
        yield (pred, _get_edge_sign(im, pred_edge) * polarity)


def _get_edge_sign(im, edge):
    """Get the polarity of the influence by examining the edge sign."""
    edge_data = im[edge[0]][edge[1]]
    # Handle possible multiple edges between nodes
    signs = list(set([v['sign'] for v in edge_data.values()
                                  if v.get('sign')]))
    if len(signs) > 1:
        logger.warning("Edge %s has conflicting polarities; choosing "
                       "positive polarity by default" % str(edge))
        sign = 1
    else:
        sign = signs[0]
    if sign is None:
        raise Exception('No sign attribute for edge.')
    elif abs(sign) == 1:
        return sign
    else:
        raise Exception('Unexpected edge sign: %s' % edge.attr['sign'])


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


def _add_activity_to_agent(agent, act_type, is_active):
    # Default to active, and return polarity if it's an inhibition
    new_act = ActivityCondition(act_type, True)
    # Check if this state already exists
    if agent.activity is not None and agent.activity.equals(new_act):
        return agent
    new_agent = deepcopy(agent)
    new_agent.activity = new_act
    polarity = 1 if is_active else -1
    return (new_agent, polarity)


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
    # Check that any state in cp2 is matched in cp1
    # If the thing we're matching to is just a monomer pattern, that makes
    # things easier--we just need to find the corresponding monomer pattern
    # in cp1
    if cp1 is None or cp2 is None:
        return False
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

def _flip(im, path):
    # Reverse the path and the polarities associated with each node
    rev = tuple(reversed(path))
    return _path_with_polarities(im, rev)


def _path_with_polarities(im, path):
    # This doesn't address the effect of the rules themselves on the
    # observables of interest--just the effects of the rules on each other
    edge_polarities = []
    path_list = list(path)
    edges = zip(path_list[0:-1], path_list[1:])
    for from_tup, to_tup in edges:
        from_rule = from_tup[0]
        to_rule = to_tup[0]
        edge = (from_rule, to_rule)
        edge_polarities.append(_get_edge_sign(im, edge))
    # Compute and return the overall path polarity
    #path_polarity = np.prod(edge_polarities)
    # Calculate left product of edge polarities return
    polarities_lprod = [1]
    for ep_ix, ep in enumerate(edge_polarities):
        polarities_lprod.append(polarities_lprod[-1] * ep)
    assert len(path) == len(polarities_lprod)
    return tuple(zip([node for node, sign in path], polarities_lprod))
    #assert path_polarity == 1 or path_polarity == -1
    #return True if path_polarity == 1 else False
    #return path_polarity


def stmt_from_rule(rule_name, model, stmts):
    """Return the source INDRA Statement corresponding to a rule in a model.

    Parameters
    ----------
    rule_name : str
        The name of a rule in the given PySB model.
    model : pysb.core.Model
        A PySB model which contains the given rule.
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements from which the model was assembled.

    Returns
    -------
    stmt : indra.statements.Statement
        The Statement from which the given rule in the model was obtained.
    """
    stmt_uuid = None
    for ann in model.annotations:
        if ann.subject == rule_name:
            if ann.predicate == 'from_indra_statement':
                stmt_uuid = ann.object
                break
    if stmt_uuid:
        for stmt in stmts:
            if stmt.uuid == stmt_uuid:
                return stmt


def _monomer_pattern_label(mp):
    """Return a string label for a MonomerPattern."""
    site_strs = []
    for site, cond in mp.site_conditions.items():
        if isinstance(cond, tuple) or isinstance(cond, list):
            assert len(cond) == 2
            if cond[1] == WILD:
                site_str = '%s_%s' % (site, cond[0])
            else:
                site_str = '%s_%s%s' % (site, cond[0], cond[1])
        elif isinstance(cond, numbers.Real):
            continue
        else:
            site_str = '%s_%s' % (site, cond)
        site_strs.append(site_str)
    return '%s_%s' % (mp.monomer.name, '_'.join(site_strs))


def _im_to_signed_digraph(im):
    edges = []
    for e in im.edges():
        edge_sign = _get_edge_sign(im, e)
        polarity = 0 if edge_sign > 0 else 1
        edges.append((e[0], e[1], dict([('sign', polarity)])))
    dg = nx.DiGraph()
    dg.add_edges_from(edges)
    return dg

