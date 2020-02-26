import logging
import textwrap
import itertools
from collections import deque
from copy import deepcopy

import numpy as np
import networkx as nx

try:
    import paths_graph as pg
    has_pg = True
except ImportError:
    has_pg = False


logger = logging.getLogger(__name__)


class PathMetric(object):
    """Describes results of simple path search (path existence).

    Attributes
    ----------
    source_node : str
        The source node of the path
    target_node : str
        The target node of the path
    length : int
        The length of the path
    """
    def __init__(self, source_node, target_node, length):
        self.source_node = source_node
        self.target_node = target_node
        self.length = length

    def __repr__(self):
        return str(self)

    def __str__(self):
        return ('source_node: %s, target_node: %s, length: %d' %
                (self.source_node, self.target_node, self.length))


class PathResult(object):
    """Describes results of running the ModelChecker on a single Statement.

    Attributes
    ----------
    path_found : bool
        True if a path was found, False otherwise.
    result_code : string
        - *STATEMENT_TYPE_NOT_HANDLED* - The provided statement type is not
          handled
        - *SUBJECT_MONOMERS_NOT_FOUND* or *SUBJECT_NOT_FOUND* -
          Statement subject not found in model
        - *OBSERVABLES_NOT_FOUND* or *OBJECT_NOT_FOUND* -
          Statement has no associated observable
        - *NO_PATHS_FOUND* - Statement has no path for any observable
        - *MAX_PATH_LENGTH_EXCEEDED* - Statement has no path len <=
          MAX_PATH_LENGTH
        - *PATHS_FOUND* - Statement has path len <= MAX_PATH_LENGTH
        - *INPUT_RULES_NOT_FOUND* - No rules with Statement subject found
        - *MAX_PATHS_ZERO* - Path found but MAX_PATHS is set to zero
    max_paths : int
        The maximum number of specific paths to return for each Statement
        to be explained.
    max_path_length : int
        The maximum length of specific paths to return.
    path_metrics : list[:py:class:`indra.explanation.model_checker.PathMetric`]
        A list of PathMetric objects, each describing the results of a simple
        path search (path existence).
    paths : list[list[tuple[str, int]]]
        A list of paths obtained from path finding. Each path is a list of
        tuples (which are edges in the path), with the first element of the
        tuple the name of a rule, and the second element its polarity in the
        path.
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
    """The parent class of all ModelCheckers.

    Parameters
    ----------
    model : pysb.Model or indra.assemblers.indranet.IndraNet or PyBEL.Model
        Depending on the ModelChecker class, can be different type.
    statements : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to check the model against.
    do_sampling : bool
        Whether to use breadth-first search or weighted sampling to
        generate paths. Default is False (breadth-first search).
    seed : int
        Random seed for sampling (optional, default is None).

    Attributes
    ----------
    graph : nx.Digraph
        A DiGraph with signed nodes to find paths in.
    """
    def __init__(self, model, statements=None, do_sampling=False, seed=None):
        self.model = model
        if statements:
            self.statements = statements
        else:
            self.statements = []
        if seed is not None:
            np.random.seed(seed)
        # Whether to do sampling
        self.do_sampling = do_sampling
        self.graph = None

    def add_statements(self, stmts):
        """Add to the list of statements to check against the model.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            The list of Statements to be added for checking.
        """
        self.statements += stmts

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
        for idx, stmt in enumerate(self.statements):
            logger.info('---')
            logger.info('Checking statement (%d/%d): %s' %
                        (idx + 1, len(self.statements), stmt))
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
        result : indra.explanation.modelchecker.PathResult
            A PathResult object containing the result of a test.
        """
        # Make sure graph is created
        self.get_graph()
        # Extract subject and object info from test statement
        subj_list, obj_list, result_code = self.process_statement(stmt)
        if result_code:
            return self.make_false_result(result_code, max_paths,
                                          max_path_length)
        # This is the case if we are checking a Statement whose
        # subject is genuinely None
        if all(s is None for s in subj_list):
            input_set = None
        # This is the case where the Statement has an actual subject
        # but we may still run into issues with finding an input
        # set for it in which case a false result may be returned.
        else:
            logger.info('Subject list: %s' % str(subj_list))
            input_set = []
            meaningful_res_code = None
            # Each subject might produce a different input set and we need to
            # combine them
            for subj in subj_list:
                inp, res_code = self.process_subject(subj)
                if res_code:
                    meaningful_res_code = res_code
                    continue
                input_set += inp
            if not input_set and meaningful_res_code:
                return self.make_false_result(meaningful_res_code, max_paths,
                                              max_path_length)

        logger.info('Input set: %s' % str(input_set))

        # If source and target are the same, we need to handle a loop
        loop = False
        if (input_set and (len(input_set) == len(obj_list) == 1) and
                (list(input_set)[0] == list(obj_list)[0])):
            loop = True
        # Now we add a dummy target node as a child to all nodes in obj_list
        common_target = ('common_target', 0)
        self.graph.add_node(common_target)
        # This is the case when source and target are the same. NetworkX does
        # not allow loops in the paths, so we work around it by using target
        # predecessors as new targets
        if loop:
            for obj in self.graph.predecessors(list(obj_list)[0]):
                self.graph.add_edge(obj, common_target)
        else:
            for obj in obj_list:
                self.graph.add_edge(obj, common_target)
        result = self.find_paths(input_set, common_target, max_paths,
                                 max_path_length, loop)

        self.graph.remove_node(common_target)
        # If a path was found, then we return it; otherwise, that means
        # there was no path for this object, so we have to try the next
        # one
        if result.path_found:
            logger.info('Found paths for %s' % stmt)
            return result
        # Return the result if the subject/input rules were not found
        if result.result_code in [
                'SUBJECT_NOT_FOUND', 'INPUT_RULES_NOT_FOUND']:
            return result
        # If we got here, then there was no path for any observable
        logger.info('No paths found for %s' % stmt)
        return self.make_false_result('NO_PATHS_FOUND',
                                      max_paths, max_path_length)

    def find_paths(self, input_set, target, max_paths=1, max_path_length=5,
                   loop=False):
        """Check for a source/target path in the model.

        Parameters
        ----------
        input_set : list or None
            A list of potenital sources or None if the test statement subject
            is None.
        target : tuple
            Tuple representing the target node (usually common target node).
        max_paths : int
            The maximum number of specific paths to return.
        max_path_length : int
            The maximum length of specific paths to return.
        loop : bool
            Whether we are looking for a loop path.

        Returns
        -------
        PathResult
            PathResult object indicating the results of the attempt to find
            a path.
        """
        # # -- Route to the path sampling function --
        # NOTE this is not generic at this point!
        # if self.do_sampling:
        #     if not has_pg:
        #         raise Exception('The paths_graph package could not be '
        #                         'imported.')
        #     return self._sample_paths(input_set, obj, target_polarity,
        #                               max_paths, max_path_length)

        # -- Do Breadth-First Enumeration --
        # Generate the predecessors to our observable and count the paths
        path_lengths = []
        path_metrics = []
        sources = []
        for source, path_length in self._find_sources(target, input_set):
            # Path already includes an edge from targets to common target, so
            # we need to subtract one edge. In case of loops, we are
            # already missing one edge, there's no need to subtract one more.
            if not loop:
                path_length = path_length - 1
            # There might be a case when sources and targets contain the same
            # nodes (e.g. different agent state in PyBEL networks) that would
            # show up as paths of length 0. We only want to include meaningful
            # paths that contain at least one edge.
            if path_length > 0:
                pm = PathMetric(source, target, path_length)
                path_metrics.append(pm)
                path_lengths.append(path_length)
                # Keep unique sources but use a list (not set) to preserve order
                if source not in sources:
                    sources.append(source)
        # Now, look for paths
        if path_metrics and max_paths == 0:
            pr = PathResult(True, 'MAX_PATHS_ZERO',
                            max_paths, max_path_length)
            pr.path_metrics = path_metrics
            return pr
        elif path_metrics:
            if min(path_lengths) <= max_path_length:
                if not loop:
                    search_path_length = min(path_lengths) + 1
                else:
                    search_path_length = min(path_lengths)
                pr = PathResult(True, 'PATHS_FOUND',
                                max_paths, max_path_length)
                pr.path_metrics = path_metrics
                # Get the first path
                # Try to find paths of fixed length using sources found above
                for source in sources:
                    logger.info('Finding paths between %s and %s'
                                % (str(source), target))
                    path_iter = get_path_iter(self.graph, source, target,
                                              search_path_length, loop)
                    for path in path_iter:
                        pr.add_path(tuple(path))
                        # Do not get next path if reached max_paths
                        if len(pr.paths) >= max_paths:
                            break
                    # Do not check next source if reached max_paths
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

    def _find_sources(self, target, sources):
        """Get the set of source nodes with paths to the target.

        Given a common target and  a list of sources (or None if test statement
        subject is None), perform a breadth-first search upstream from the
        target to determine whether there are any sources that have paths to
        the target. For efficiency, does not return the full path,
        but identifies the upstream sources and the length of the path.

        Parameters
        ----------
        target : tuple
            The node (usually common target node) in the graph to start
            looking upstream for matching sources.
        sources : list[tuple]
            Signed nodes corresponding to the subject or upstream influence
            being checked.

        Returns
        -------
        generator of (source, path_length)
            Yields tuples of source node (tuple of string and sign) and path
            length (int). If there are no paths to any of the given source
            nodes, the generator is empty.
        """
        # First, create a list of visited nodes
        # Adapted from
        # networkx.algorithms.traversal.breadth_first_search.bfs_edges
        visited = set([target])
        # Generate list of predecessor nodes with a sign updated according to
        # the sign of the target node

        # The queue holds tuples of "parents" (in this case downstream nodes)
        # and their "children" (in this case their upstream influencers)
        queue = deque([(target, self.graph.predecessors(target), 0)])
        while queue:
            parent, children, path_length = queue[0]
            try:
                # Get the next child in the list
                child = next(children)
                # Is this child one of the source nodes we're looking for? If
                # so, yield it along with path length.
                # Also make sure that found source is positive
                if (sources is None or child in sources) and child[1] == 0:
                    logger.debug("Found path to %s from %s with length %d"
                                 % (target, child, path_length+1))
                    yield (child, path_length+1)
                # Check this child against the visited list. If we haven't
                # visited it already (accounting for the path to the node),
                # then add it to the queue.
                if child not in visited:
                    visited.add(child)
                    queue.append(
                        (child, self.graph.predecessors(child), path_length + 1))
            # Once we've finished iterating over the children of the current
            # node, pop the node off and go to the next one in the queue
            except StopIteration:
                queue.popleft()
        # There was no path; this will produce an empty generator
        return

    def make_false_result(self, result_code, max_paths, max_path_length):
        return PathResult(False, result_code, max_paths, max_path_length)

    def get_graph(self, **kwargs):
        """Return a graph  with signed nodes to find the path."""
        raise NotImplementedError("Method must be implemented in child class.")

    def process_statement(self, stmt):
        """
        This method processes the test statement to get the data about subject
        and object, according to the specific model requirements for model
        checking, e.g. PysbModelChecker gets subject monomer patterns and
        observables, while graph based ModelCheckers will return signed nodes
        corresponding to subject and object. If any of the requirements are not
        satisfied, result code is also returned to construct PathResult object.

        Parameters
        ----------
        stmt : indra.statements.Statement
            A statement to process.

        Returns
        -------
        subj_data : list or None
            Data about statement subject to be used as source nodes.
        obj_data : list or None
            Data about statement object to be used as target nodes.
        result_code : str or None
            Result code to construct PathResult.
        """
        raise NotImplementedError("Method must be implemented in child class.")

    def process_subject(self, subj_data):
        """Processes the subject of the test statement and returns
        the necessary information to check the statement. In case of
        PysbModelChecker, method returns input_rule_set. If any of the
        requirements are not satisfied, result code is also returned to
        construct PathResult object.
        """
        raise NotImplementedError("Method must be implemented in child class.")

    def _sample_paths(self, input_set, obj_name, target_polarity,
                      max_paths=1, max_path_length=5):
        raise NotImplementedError("Method must be implemented in child class.")


def get_path_iter(graph, source, target, path_length, loop):
    """Return a generator of paths with path_length cutoff from source to target."""
    path_iter = nx.all_simple_paths(graph, source, target, path_length)
    try:
        for p in path_iter:
            path = deepcopy(p)
            # Remove common target from a path.
            path.remove(target)
            if loop:
                path.append(path[0])
            # A path should contain at least one edge
            if len(path) < 2:
                continue
            yield path
    except nx.NetworkXNoPath:
        pass


def signed_edges_to_signed_nodes(graph, prune_nodes=True,
                                 edge_signs={'pos': 0, 'neg': 1},
                                 copy_edge_data=False):
    """Convert a graph with signed edges to a graph with signed nodes.

    Each pair of nodes linked by an edge in an input graph are represented
    as four nodes and two edges in the new graph. For example, an edge (a,
    b, 0), where a and b are nodes and 0 is a sign of an edge (positive),
    will be represented as edges ((a, 0), (b, 0)) and ((a, 1), (b, 1)),
    where (a, 0), (a, 1), (b, 0), (b, 1) are signed nodes. An edge (a, b,
    1) with sign 1 (negative) will be represented as edges ((a, 0), (b,
    1)) and ((a, 1), (b, 0)).

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        Graph with signed edges to convert. Can have multiple edges between
        a pair of nodes.
    prune_nodes : Optional[bool]
        If True, iteratively prunes negative (with sign 1) nodes without
        predecessors.
    edge_signs : dict
        A dictionary representing the signing policy of incoming graph. The
        dictionary should have strings 'pos' and 'neg' as keys and integers
        as values.
    copy_edge_data : bool|set(keys)
        Option for copying edge data as well from graph. If False (default),
        no edge data is copied (except sign). If True, all edge data is
        copied. If a set of keys is provided, only the keys appearing in the
        set will be copied, assuming the key is part of a nested dictionary.

    Returns
    -------
    signed_nodes_graph : networkx.DiGraph
    """
    signed_nodes_graph = nx.DiGraph()
    nodes = []
    for node, node_data in graph.nodes(data=True):
        nodes.append(((node, 0), node_data))
        nodes.append(((node, 1), node_data))
    signed_nodes_graph.add_nodes_from(nodes)
    edges = []
    for u, v, edge_data in graph.edges(data=True):
        copy_dict = deepcopy(edge_data)
        edge_sign = copy_dict.pop('sign', None)
        if edge_sign is None:
            continue
        edge_dict = copy_dict if copy_edge_data == True else \
            ({k: v for k, v in copy_dict.items() if k in copy_edge_data} if
             isinstance(copy_edge_data, set) else {})
        if edge_sign == edge_signs['pos']:
            edges.append(((u, 0), (v, 0), edge_dict))
            edges.append(((u, 1), (v, 1), edge_dict))
        elif edge_sign == edge_signs['neg']:
            edges.append(((u, 0), (v, 1), edge_dict))
            edges.append(((u, 1), (v, 0), edge_dict))
    signed_nodes_graph.add_edges_from(edges)
    if prune_nodes:
        signed_nodes_graph = prune_signed_nodes(signed_nodes_graph)
    return signed_nodes_graph


def prune_signed_nodes(graph):
    """Prune nodes with sign (1) if they do not have predecessors."""
    nodes_to_prune = [node for node, in_deg
                      in graph.in_degree()
                      if in_deg == 0 and node[1] == 1]
    while nodes_to_prune:
        graph.remove_nodes_from(nodes_to_prune)
        # Make a list of nodes whose in degree is now 0
        nodes_to_prune = [node for node, in_deg
                          in graph.in_degree()
                          if in_deg == 0 and node[1] == 1]
    return graph
