__all__ = ['get_agent_key', 'get_relevant_keys', 'RefinementFilter',
           'RefinementConfirmationFilter', 'OntologyRefinementFilter',
           'ontology_refinement_filter', 'bio_ontology_refinement_filter',
           'default_refinement_fun']

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
    def __init__(self):
        self.shared_data = {}

    def initialize(self, stmts_by_hash):
        self.shared_data['stmts_by_hash'] = stmts_by_hash

    def get_more_specifics(self, stmt, possibly_related=None):
        pass

    def get_less_specifics(self, stmt, possibly_related=None):
        pass


class OntologyRefinementFilter(RefinementFilter):
    def __init__(self, ontology):
        self.ontology = ontology
        self.shared_data = {}

    def initialize(self, stmts_by_hash):
        self.shared_data['stmts_by_hash'] = stmts_by_hash

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
            roles = stmts_by_hash[
                next(iter(stmts_this_type))]._agent_order
            # Mapping agent keys to statement hashes
            agent_key_to_hash = {}
            # Mapping statement hashes to agent keys
            hash_to_agent_key = {}
            # All agent keys for a given agent role
            all_keys_by_role = {}
            for role in roles:
                agent_key_to_hash[role] = collections.defaultdict(set)
                hash_to_agent_key[role] = collections.defaultdict(set)

            # Step 2. Fill up the initial data structures in preparation
            # for identifying potential refinements
            for sh in stmts_this_type:
                for role in roles:
                    agent_keys = self._agent_keys_for_stmt_role(
                        stmts_by_hash[sh], role)
                    for agent_key in agent_keys:
                        agent_key_to_hash[role][agent_key].add(sh)
                        hash_to_agent_key[role][sh].add(agent_key)

            agent_key_to_hash = dict(agent_key_to_hash)
            hash_to_agent_key = dict(hash_to_agent_key)

            for role in roles:
                all_keys_by_role[role] = set(agent_key_to_hash[role].keys())

            # Step 3. Make these available as shared data for the apply step
            self.shared_data[stmt_type] = {
                'agent_key_to_hash': agent_key_to_hash,
                'hash_to_agent_key': hash_to_agent_key,
                'all_keys_by_role': all_keys_by_role,
            }

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

    def get_less_specifics(self, stmt, possibly_related=None):
        return self._get_related(stmt, possibly_related=possibly_related,
                                 direction='less_specific')

    def get_more_specifics(self, stmt, possibly_related=None):
        return self._get_related(stmt, possibly_related=possibly_related,
                                 direction='more_specific')

    def _get_related(self, stmt, possibly_related=None,
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
    def __init__(self, ontology, refinement_fun=None):
        self.ontology = ontology
        self.refinement_fun = refinement_fun if refinement_fun else \
            default_refinement_fun
        self.shared_data = {}
        self.comparison_counter = 0

    def get_less_specifics(self, stmt, possibly_related=None):
        return self._get_related(stmt, possibly_related=possibly_related,
                                 direction='less_specific')

    def get_more_specifics(self, stmt, possibly_related=None):
        return self._get_related(stmt, possibly_related=possibly_related,
                                 direction='more_specific')

    def _get_related(self, stmt, possibly_related=None,
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


def ontology_refinement_filter(stmts_by_hash, stmts_to_compare, ontology):
    """Return possible refinement relationships based on an ontology.

    Parameters
    ----------
    stmts_by_hash : dict
        A dict whose keys are statement hashes that point to the
        (deduplicated) statement with that hash as a value.
    stmts_to_compare : dict or None
        A dict of existing statements to compare that will be further
        filtered down in this function and then returned.
    ontology : indra.ontology.IndraOntology
        An IndraOntology instance iwth respect to which this
        filter is applied.

    Returns
    -------
    dict
        A dict whose keys are statement hashes and values are sets
        of statement hashes that can potentially be refined by the
        statement identified by the key.
    """
    logger.info('Finding ontology-based refinements for %d statements'
                % len(stmts_by_hash))
    ts = time.time()
    ont_filter = OntologyRefinementFilter(ontology)
    ont_filter.initialize(stmts_by_hash)

    stmts_to_compare = {
        stmt_hash: ont_filter.get_less_specifics(
            stmt_hash,
            possibly_related=(
                stmts_to_compare.get(stmt_hash)
                if stmts_to_compare is not None else None
                )
            )
        for stmt_hash in stmts_by_hash}
    te = time.time()
    logger.debug('Identified ontology-based possible refinements in %.2fs'
                 % (te-ts))
    # Make an empty dict to make sure we don't return a None
    if stmts_to_compare is None:
        stmts_to_compare = {}
    return stmts_to_compare


def bio_ontology_refinement_filter(stmts_by_hash, stmts_to_compare):
    """An ontology refinement filter that works with the INDRA BioOntology."""
    from indra.ontology.bio import bio_ontology
    return ontology_refinement_filter(stmts_by_hash, stmts_to_compare,
                                      ontology=bio_ontology)


def default_refinement_fun(st1, st2, ontology, entities_refined):
    return st1.refinement_of(st2, ontology, entities_refined)
