"""This module implements classes and functions that are used for
finding refinements between INDRA Statements as part of the
knowledge-assembly process. These are imported by the preassembler
module."""
__all__ = ['get_agent_key', 'get_relevant_keys', 'RefinementFilter',
           'RefinementConfirmationFilter', 'OntologyRefinementFilter',
           'SplitGroupFilter', 'default_refinement_fun']

import time
import logging
import collections
from indra.statements import Event
from indra.statements import stmt_type as indra_stmt_type


logger = logging.getLogger(__name__)


# TODO: we could make the agent key function parameterizable with the
# preassembler to allow custom agent mappings to the ontology.
def get_agent_key(agent):
    """Return a key for an Agent for use in refinement finding.

    Parameters
    ----------
    agent : indra.statements.Agent or None
         An INDRA Agent whose key should be returned.

    Returns
    -------
    tuple or None
        The key that maps the given agent to the ontology, with special
        handling for ungrounded and None Agents.
    """
    if isinstance(agent, Event):
        agent = agent.concept
    if agent is None:
        agent_key = None
    else:
        agent_key = agent.get_grounding()
        if not agent_key[0]:
            agent_key = ('NAME', agent.name)
    return agent_key


def get_relevant_keys(agent_key, all_keys_for_role, ontology, direction):
    """Return relevant agent keys for an agent key for refinement finding.

    Parameters
    ----------
    agent_key : tuple or None
        An agent key of interest.
    all_keys_for_role : set
        The set of all agent keys in a given statement corpus with a
        role matching that of the given agent_key.
    ontology : indra.ontology.IndraOntology
        An IndraOntology instance with respect to which relevant other
        agent keys are found for the purposes of refinement.
    direction: str
        The direction in which to find relevant agents. The two options
        are 'less_specific' and 'more_specific' for agents that are less and
        more specific, per the ontology, respectively.

    Returns
    -------
    set
        The set of relevant agent keys which this given agent key can
        possibly refine.
    """
    rel_fun = ontology.get_parents if direction == 'less_specific' else \
        ontology.get_children
    relevant_keys = {None, agent_key}
    if agent_key is not None:
        relevant_keys |= set(rel_fun(*agent_key))
    relevant_keys &= all_keys_for_role
    return relevant_keys


class RefinementFilter:
    """A filter which is applied to one or more statements to eliminate
    candidate refinements that are not possible according to some
    criteria. By applying a series of such filters, the preassembler can avoid
    doing n-by-n comparisons to determine refinements among n statements.

    The filter class can take any number of constructor arguments that it
    needs to perform its task. The base class' constructor initializes
    a shared_data attribute as an empty dict.

    It also needs to implement an initialize function which is called
    with a stmts_by_hash argument, containing a dict of statements keyed by
    hash. This function can build any data structures that may be needed
    to efficiently apply the filter later. It cab store any
    such data structures in the shared_data dict to be accessed by
    other functions later.

    Finally, the class needs to implement a get_related function, which
    takes a single INDRA Statement as input to return the hashes of
    potentially related other statements that the filter was initialized
    with. The function also needs to take a possibly_related argument
    which is either None (no other filter was run before) or a set,
    which is the superset of possible relations as determined by some
    other previously applied filter.
    """
    def __init__(self):
        self.shared_data = {}

    def initialize(self, stmts_by_hash):
        """Initialize the filter class with a set of statements.

        The filter can build up some useful data structures in this
        function before being applied to any specific statements.

        Parameters
        ----------
        stmts_by_hash : dict[int, indra.statements.Statement]
            A dict of statements keyed by their hashes.
        """
        self.shared_data['stmts_by_hash'] = stmts_by_hash

    def get_related(self, stmt, possibly_related=None,
                    direction='less_specific'):
        """Return a set of statement hashes that a given statement is
        potentially related to.

        Parameters
        ----------
        stmt : indra.statements.Statement
            The INDRA statement whose potential relations we want to filter.
        possibly_related : set or None
            A set of statement hashes that this statement is potentially
            related to, as determined by some other filter. If this parameter
            is a set (including an empty set), this function should return
            a subset of it (intuitively, this filter can only further eliminate
            some of the potentially related hashes that were previously
            determined to be potential relations). If this argument is
            None, the function must assume that no previous filter
            was run before, and should therefore return all the possible
            relations that it determines.
        direction : str
            One of 'less_specific' or 'more_specific. Since refinements
            are directed relations, this function can operate in two
            different directions: it can either find less specific
            potentially related stateemnts, or it can find more specific
            potentially related statements, as determined by this argument.

        Returns
        -------
        set of int
            A set of INDRA Statement hashes that are potentially related
            to the given statement.
        """
        raise NotImplementedError('The filter class has to implement a'
                                  'get_related method.')

    def get_more_specifics(self, stmt, possibly_related=None):
        """Return a set of hashes of statements that are potentially related
        and more specific than the given statement."""
        return self.get_related(stmt, possibly_related=possibly_related,
                                direction='more_specific')

    def get_less_specifics(self, stmt, possibly_related=None):
        """Return a set of hashes of statements that are potentially related
        and less specific than the given statement."""
        return self.get_related(stmt, possibly_related=possibly_related,
                                direction='less_specific')

    def extend(self, stmts_by_hash):
        """Extend the initial data structures with a set of new statements.

        Parameters
        ----------
        stmts_by_hash : dict[int, indra.statements.Statement]
            A dict of statements keyed by their hashes.
        """
        # We can assume that these stmts_by_hash are unique
        self.shared_data['stmts_by_hash'].update(stmts_by_hash)


class OntologyRefinementFilter(RefinementFilter):
    """This filter uses an ontology to position statements and their agents
    to filter down significantly on the set of possible relations for
    a given statement.

    Parameters
    ----------
    ontology : indra.ontology.OntologyGraph
        An INDRA ontology graph.
    """
    def __init__(self, ontology):
        super().__init__()
        self.ontology = ontology

    def initialize(self, stmts_by_hash):
        self.shared_data['stmts_by_hash'] = {}
        self.extend(stmts_by_hash)

    def extend(self, stmts_by_hash):
        self.shared_data['stmts_by_hash'].update(stmts_by_hash)
        # Build up data structure of statement hashes by
        # statement type
        stmts_by_type = collections.defaultdict(set)
        for stmt_hash, stmt in stmts_by_hash.items():
            stmts_by_type[indra_stmt_type(stmt)].add(stmt_hash)
        stmts_by_type = dict(stmts_by_type)

        # Now iterate over each statement type and build up
        # data structures for quick filtering
        for stmt_type, stmts_this_type in stmts_by_type.items():
            # Step 1. initialize data structures
            # noinspection PyProtectedMember
            roles = stmts_by_hash[next(iter(stmts_this_type))]._agent_order
            if stmt_type not in self.shared_data:
                self.shared_data[stmt_type] = {}
                # Mapping agent keys to statement hashes
                self.shared_data[stmt_type]['agent_key_to_hash'] = \
                    {role: collections.defaultdict(set) for role in roles}
                # Mapping statement hashes to agent keys
                self.shared_data[stmt_type]['hash_to_agent_key'] = \
                    {role: collections.defaultdict(set) for role in roles}
                # All agent keys for a given agent role
                self.shared_data[stmt_type]['all_keys_by_role'] = {}

            # Step 2. Fill up the initial data structures in preparation
            # for identifying potential refinements
            for sh in stmts_this_type:
                for role in roles:
                    agent_keys = self._agent_keys_for_stmt_role(
                        stmts_by_hash[sh], role)
                    for agent_key in agent_keys:
                        self.shared_data[stmt_type]['agent_key_to_hash'][
                            role][agent_key].add(sh)
                        self.shared_data[stmt_type]['hash_to_agent_key'][
                            role][sh].add(agent_key)

            for role in roles:
                self.shared_data[stmt_type]['all_keys_by_role'][role] = \
                    set(self.shared_data[stmt_type]['agent_key_to_hash'][role])

    @staticmethod
    def _agent_keys_for_stmt_role(stmt, role):
        """Return a set of agent keys for a statement's agent in a role.

        The agent key is an "anchor" to the ontology being used and positons
        a statement, via its agent in this role against other statements it
        may be related to.
        """
        agents = getattr(stmt, role)
        # Handle a special case here where a list=like agent
        # role can be empty, here we will consider anything else
        # to be a refinement, hence add a None key
        if isinstance(agents, list) and not agents:
            agent_keys = {None}
        # Generally, we take all the agent keys for a single or
        # list-like agent role.
        else:
            agent_keys = {get_agent_key(agent) for agent in
                          (agents if isinstance(agents, list)
                           else [agents])}
        return agent_keys

    def get_related(self, stmt, possibly_related=None,
                    direction='less_specific'):
        # Corner case: if this is a new statement that wasn't part of the
        # initialization, it is possible that it has a type that we've not
        # seen during initialization at all. In this case, we can assume
        # there are no refinements for it.
        stmt_type = indra_stmt_type(stmt)
        if stmt_type not in self.shared_data:
            return {}

        # Step 1. Recover relevant parts ot the initialized data
        hash_to_agent_key = self.shared_data[stmt_type]['hash_to_agent_key']
        agent_key_to_hash = self.shared_data[stmt_type]['agent_key_to_hash']
        all_keys_by_role = self.shared_data[stmt_type]['all_keys_by_role']

        # Step 2. We iterate over all statements and find ones that this one
        # can refine
        stmt_hash = stmt.get_hash()
        relevants = possibly_related
        # We now iterate over all the agent roles in the given statement
        # type
        for role, hash_to_agent_key_for_role in hash_to_agent_key.items():
            # If we have seen this statement before during initialization then
            # we can use its precalculated agent keys, otherwise we
            # calculate new agent keys for it.
            if stmt_hash in hash_to_agent_key_for_role:
                agent_keys = hash_to_agent_key_for_role[stmt_hash]
            else:
                agent_keys = self._agent_keys_for_stmt_role(stmt, role)

            # We get all the agent keys in all other statements that the
            # agent in this given role in this statement can refine.
            for agent_key in agent_keys:
                relevant_keys = get_relevant_keys(
                    agent_key,
                    all_keys_by_role[role],
                    self.ontology,
                    direction=direction)
                # We now get the actual statement hashes that these other
                # potentially refined agent keys appear in in the given role
                role_relevant_stmt_hashes = set.union(
                    *[agent_key_to_hash[role][rel]
                      for rel in relevant_keys]) - {stmt_hash}
                # In the first iteration, we initialize the set with the
                # relevant statement hashes
                if relevants is None:
                    relevants = role_relevant_stmt_hashes
                # In subsequent iterations, we take the intersection of
                # the relevant sets per role
                else:
                    relevants &= role_relevant_stmt_hashes

        # These hashes are now the ones that this statement needs
        # to be compared against. Importantly, the relationship is in
        # a well-defined direction so we don't need to test both ways.
        return relevants


class RefinementConfirmationFilter(RefinementFilter):
    """This class runs the refinement function between potentially
    related statements to confirm whether they are indeed, conclusively
    in a refinement relationship with each other.

    In this sense, this isn't a real filter, though implementing it
    as one is convenient. This filter is meant to be used as the final
    component in a series of pre-filters.
    """
    def __init__(self, ontology, refinement_fun=None):
        self.ontology = ontology
        self.refinement_fun = refinement_fun if refinement_fun else \
            default_refinement_fun
        self.shared_data = {}
        self.comparison_counter = 0

    def get_related(self, stmt, possibly_related=None,
                    direction='less_specific'):
        stmts_by_hash = self.shared_data['stmts_by_hash']
        relateds = set()
        # We use the previously constructed set of statements that this one
        # can possibly refine
        for possible_related_hash in possibly_related:
            more_spec_stmt, less_spec_stmt = (
                (stmt, stmts_by_hash[possible_related_hash])
                if direction == 'less_specific'
                else (stmts_by_hash[possible_related_hash], stmt)
            )
            # And then do the actual comparison. Here we use
            # entities_refined=True which means that we assert that
            # the entities, in each role, are already confirmed to
            # be "compatible" for refinement, and therefore, we
            # don't need to again confirm this (i.e., call "isa") in
            # the refinement_of function.
            ref = self.refinement_fun(
                more_spec_stmt,
                less_spec_stmt,
                ontology=self.ontology,
                # NOTE: here we assume that the entities at this point
                # are definitely refined due to the use of an
                # ontology-based pre-filter. If this is not the case
                # for some reason then it is the responsibility of the
                # user-supplied self.refinement_fun to disregard the
                # entities_refined argument.
                entities_refined=True)
            self.comparison_counter += 1
            if ref:
                relateds.add(possible_related_hash)
        return relateds


class SplitGroupFilter(RefinementFilter):
    """This filter implements splitting statements into two groups and
    only considering refinement relationships between the groups but not
    within them."""
    def __init__(self, split_groups):
        super().__init__()
        self.split_groups = split_groups

    def get_related(self, stmt, possibly_related=None,
                    direction='less_specific'):
        sh = stmt.get_hash()
        group = self.split_groups.get(sh)
        # We take all the hashes that are in a different group, and
        # return all of them of possibly_related is None (i.e., there
        # was no previous filter), or if the given statement is
        # also possibly related
        related = {stmt_hash for stmt_hash
                   in self.shared_data['stmts_by_hash']
                   if self.split_groups[stmt_hash] != group
                   and (possibly_related is None
                        or stmt_hash in possibly_related)}
        return related


def default_refinement_fun(st1, st2, ontology, entities_refined):
    return st1.refinement_of(st2, ontology, entities_refined)
