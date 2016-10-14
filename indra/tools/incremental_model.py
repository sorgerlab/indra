from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import logging
from indra.statements import Agent
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler import grounding_mapper as gm
from indra.belief import BeliefEngine
from indra.databases import uniprot_client

logger = logging.getLogger('incremental_model')

class IncrementalModel(object):
    """Assemble a model incrementally by iteratively adding new Statements.

    Parameters
    ----------
    model_fname : Optional[str]
        The name of the pickle file in which a set of INDRA Statements are
        stored in a dict keyed by PubMed IDs. This is the state of an
        IncrementalModel that is loaded upon instantiation.

    Attributes
    ----------
    stmts : dict[str, list[indra.statements.Statement]]
        A dictionary of INDRA Statements keyed by PMIDs that stores the current
        state of the IncrementalModel.
    unique_stmts : list[indra.statements.Statement]
        A list of INDRA Statements after de-duplication.
    toplevel_stmts : list[indra.statements.Statement]
        A list of the top level (most specific) INDRA Statements after
        preassembly.
    """
    def __init__(self, model_fname=None):
        if model_fname is None:
            self.stmts = {}
        else:
            try:
                with open(model_fname, 'rb') as f:
                    self.stmts = pickle.load(f)
            except:
                logger.warning('Could not load %s, starting new model.' %
                               model_fname)
                self.stmts = {}
        self.relevant_stmts = []
        self.unique_stmts = []
        self.toplevel_stmts = []

    def save(self, model_fname='model.pkl'):
        """Save the state of the IncrementalModel in a pickle file.

        model_fname : Optional[str]
            The name of the pickle file to save the state of the
            IncrementalModel in. Default: model.pkl
        """
        with open(model_fname, 'wb') as fh:
            pickle.dump(self.stmts, fh, protocol=2)

    def _relevance_filter(self, stmts, filters=None):
        stmts_to_add = list(range(len(stmts)))
        # Filter for grounding
        if 'grounding' in filters:
            for i, stmt in enumerate(stmts):
                # Check that all agents are grounded
                for ag in stmt.agent_list():
                    if ag is None:
                        continue
                    if not ag.db_refs or ag.db_refs.keys() == ['TEXT']:
                        stmts_to_add.remove(i)
                        break

        if ('prior_all' in filters) or ('prior_one' in filters):
            prior_agents = self.get_prior_agents()
            if prior_agents:
                for i in stmts_to_add:
                    agents = [a for a in stmts[i].agent_list()
                              if a is not None]
                    if 'prior_all' in filters:
                        for st_agent in agents:
                            found = False
                            for pr_agent in prior_agents:
                                if self._agent_related(st_agent, pr_agent):
                                    found = True
                                    break
                            if not found:
                                stmts_to_add.remove(i)
                                break
                    if 'prior_one' in filters:
                        found = False
                        for st_agent in agents:
                            for pr_agent in prior_agents:
                                if self._agent_related(st_agent, pr_agent):
                                    found = True
                                    break
                            if found:
                                break
                        if not found:
                                stmts_to_add.remove(i)

        if ('model_all' in filters) or ('model_one' in filters):
            model_agents = self.get_model_agents()
            if model_agents:
                for i in stmts_to_add:
                    agents = [a for a in stmts[i].agent_list()
                              if a is not None]
                    if 'model_all' in filters:
                        for st_agent in agents:
                            found = False
                            for mo_agent in model_agents:
                                if self._agent_related(st_agent, mo_agent):
                                    found = True
                                    break
                            if not found:
                                stmts_to_add.remove(i)
                                break
                    if 'model_one' in filters:
                        found = False
                        for st_agent in agents:
                            for mo_agent in model_agents:
                                if self._agent_related(st_agent, mo_agent):
                                    found = True
                                    break
                            if found:
                                break
        if ('human_only' in filters):
            for i in stmts_to_add:
                agents = [a for a in stmts[i].agent_list()
                          if a is not None]
                for a in agents:
                    up_id = a.db_refs.get('UP')
                    hgnc_id = a.db_refs.get('HGNC')
                    if not hgnc_id:
                        if up_id:
                            mnemonic = uniprot_client.get_mnemonic(up_id, True)
                            if mnemonic is None:
                                stmts_to_add.remove(i)
                                break
                            if not mnemonic.endswith('HUMAN'):
                                stmts_to_add.remove(i)
                                break

        relevant_stmts = [stmts[i] for i in stmts_to_add]
        return relevant_stmts

    def add_statements(self, pmid, stmts, filters=None):
        """Add INDRA Statements to the incremental model indexed by PMID.

        Currently the following filter options are implemented:
        - grounding: require that all Agents in statements are grounded
        - model_one: require that at least one Agent is in the incremental
                      model
        - model_all: require that all Agents are in the incremental model
        - prior_one: require that at least one Agent is in the
                      prior model
        - prior_all: require that all Agents are in the prior model
        Note that model_one -> prior_all are increasingly more restrictive
        options.

        Parameters
        ----------
        pmid : str
            The PMID of the paper from which statements were extracted.
        stmts : list[indra.statements.Statement]
            A list of INDRA Statements to be added to the model.
        filter : Optional[list[str]]
            A list of filter options to apply when adding the statements.
            See description above for more details. Default: None
        """
        # If no filter is used, we add all statements to the model
        if not filters:
            self.stmts[pmid] = stmts
            return
        # If the statements are empty in the first place
        if not stmts:
            self.stmts[pmid] = []
            return

        relevant_stmts = self._relevance_filter(stmts, filters)
        self.stmts[pmid] = relevant_stmts


    def preassemble(self, filters=None):
        """Preassemble the Statements collected in the model.

        Use INDRA's GroundingMapper, Preassembler and BeliefEngine
        on the IncrementalModel and save the unique statements and
        the top level statements in class attributes.

        Currently the following filter options are implemented:
        - grounding: require that all Agents in statements are grounded
        - model_one: require that at least one Agent is in the incremental
                      model
        - model_all: require that all Agents are in the incremental model
        - prior_one: require that at least one Agent is in the
                      prior model
        - prior_all: require that all Agents are in the prior model
        Note that model_one -> prior_all are increasingly more restrictive
        options.

        Parameters
        ----------
        filter : Optional[list[str]]
            A list of filter options to apply when choosing the statements.
            See description above for more details. Default: None
        """
        stmts = self.get_statements_noprior()

        # Fix grounding
        logger.info('Running grounding map')
        twg = gm.agent_texts_with_grounding(stmts)
        prot_map = gm.protein_map_from_twg(twg)
        gm.default_grounding_map.update(prot_map)
        gmap = gm.GroundingMapper(gm.default_grounding_map)
        gmapped_stmts = gmap.map_agents(stmts, do_rename=True)

        # Merge the prior and the mapped non-prior
        stmts = gmapped_stmts + self.get_statements_prior()
        logger.info('%d unfiltered Statements' % len(stmts))

        if filters:
            if 'grounding' in filters:
                # Filter out ungrounded statements
                logger.info('Running grounding filter')
                stmts = self._relevance_filter(stmts, ['grounding'])
                logger.info('%s Statements after filter' % len(stmts))
            if 'human_only' in filters:
                # Filter out non-human proteins
                logger.info('Running non-human protein filter')
                stmts = self._relevance_filter(stmts, ['human_only'])
                logger.info('%s Statements after filter' % len(stmts))
            for rel_key in ('prior_one', 'model_one',
                             'prior_all', 'model_all'):
                if rel_key in filters:
                    logger.info('Running relevance filter')
                    stmts = self._relevance_filter(stmts, [rel_key])
                    logger.info('%s Statements after filter' % len(stmts))

        # Combine duplicates
        logger.info('Preassembling %d Statements' % len(stmts))
        pa = Preassembler(hierarchies, stmts)
        self.unique_stmts = pa.combine_duplicates()

        # Run BeliefEngine on unique statements
        be = BeliefEngine(self.unique_stmts)
        be.set_prior_probs(self.unique_stmts)

        # Build statement hierarchy
        self.toplevel_stmts = pa.combine_related()
        # Run BeliefEngine on hierarchy
        be.set_hierarchy_probs(self.toplevel_stmts)

    def load_prior(self, prior_fname):
        """Load a set of prior statements from a pickle file.

        The prior statements have a special key in the stmts dictionary
        called "prior".
        """
        with open(prior_fname, 'rb') as fh:
            self.stmts['prior'] = pickle.load(fh)

    def get_model_agents(self):
        """Return a list of all Agents from all Statements."""
        model_stmts = self.get_statements_noprior()
        # First, get all unique name/db_refs pairs
        agents = set()
        for stmt in model_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add((a.name, a.get_grounding()))
        # Then construct sateless name/db_refs Agents
        agents = [Agent(n, db_refs={g[0]: g[1]}) for n, g in agents]
        return agents

    def get_prior_agents(self):
        """Return a list of all Agents from the prior Statements."""
        prior_stmts = self.get_statements_prior()
        # First, get all unique name/db_refs pairs
        agents = set()
        if prior_stmts is None:
            return agents
        for stmt in prior_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add((a.name, a.get_grounding()))
        agents = [Agent(n, db_refs={g[0]: g[1]}) for n, g in agents]
        return agents

    @staticmethod
    def _agent_related(a1, a2):
        if a1.matches(a2) or \
            a1.isa(a2, hierarchies) or \
            a2.isa(a1, hierarchies):
            return True
        return False

    def get_statements(self):
        """Return a list of all Statements in a single list."""
        stmt_lists = [v for k, v in self.stmts.items()]
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def get_statements_noprior(self):
        """Return a list of all non-prior Statements in a single list."""
        stmt_lists = [v for k, v in self.stmts.items() if k != 'prior']
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def get_statements_prior(self):
        """Return a list of all prior Statements in a single list."""
        if self.stmts.get('prior') is not None:
            return self.stmts['prior']
        else:
            return []
