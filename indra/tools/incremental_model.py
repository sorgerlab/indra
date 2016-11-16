from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import logging
from indra.statements import Agent
from indra.belief import BeliefEngine
from indra.assemblers import PysbAssembler
from indra.databases import uniprot_client
from indra.preassembler import Preassembler
from indra.preassembler import grounding_mapper as gm
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler.sitemapper import SiteMapper, default_site_map

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

        Parameters
        ----------
        model_fname : Optional[str]
            The name of the pickle file to save the state of the
            IncrementalModel in. Default: model.pkl
        """
        with open(model_fname, 'wb') as fh:
            pickle.dump(self.stmts, fh, protocol=2)

    def _relevance_filter(self, stmts_in, filters=None):
        if filters is None:
            filters = []
        stmts_out = stmts_in
        if 'grounding' in filters:
            stmts_out = _grounding_filter(stmts_out)
        if 'human_only' in filters:
            stmts_out = _human_only_filter(stmts_out)
        if 'model_all' in filters:
            stmts_out = _ref_agents_all_filter(stmts_out,
                                               self.get_model_agents())
        if 'model_one' in filters:
            stmts_out = _ref_agents_one_filter(stmts_out,
                                               self.get_model_agents())
        if 'prior_all' in filters:
            stmts_out = _ref_agents_all_filter(stmts_out,
                                               self.get_prior_agents())
        if 'prior_one' in filters:
            stmts_out = _ref_agents_one_filter(stmts_out,
                                               self.get_prior_agents())
        return stmts_out

    def add_statements(self, pmid, stmts, filters=None):
        """Add INDRA Statements to the incremental model indexed by PMID.

        Currently the following filter options are implemented:
        - grounding: require that all Agents in statements are grounded
        - model_one: require that at least one Agent is in the incremental model
        - model_all: require that all Agents are in the incremental model
        - prior_one: require that at least one Agent is in the prior model
        - prior_all: require that all Agents are in the prior model
        Note that model_one -> prior_all are increasingly more restrictive
        options.

        Parameters
        ----------
        pmid : str
            The PMID of the paper from which statements were extracted.
        stmts : list[indra.statements.Statement]
            A list of INDRA Statements to be added to the model.
        filters : Optional[list[str]]
            A list of filter options to apply when adding the statements.
            See description above for more details. Default: None
        """
        if pmid not in self.stmts:
            self.stmts[pmid] = []
        if not stmts:
            return
        if not filters:
            self.stmts[pmid] += stmts
            return

        relevant_stmts = self._relevance_filter(stmts, filters)
        self.stmts[pmid] += relevant_stmts

    def preassemble(self, filters=None):
        """Preassemble the Statements collected in the model.

        Use INDRA's GroundingMapper, Preassembler and BeliefEngine
        on the IncrementalModel and save the unique statements and
        the top level statements in class attributes.

        Currently the following filter options are implemented:
        - grounding: require that all Agents in statements are grounded
        - model_one: require that at least one Agent is in the incremental model
        - model_all: require that all Agents are in the incremental model
        - prior_one: require that at least one Agent is in the prior model
        - prior_all: require that all Agents are in the prior model
        Note that model_one -> prior_all are increasingly more restrictive
        options.

        Parameters
        ----------
        filters : Optional[list[str]]
            A list of filter options to apply when choosing the statements.
            See description above for more details. Default: None
        """
        stmts = self.get_statements()
        logger.info('%d raw Statements in total' % len(stmts))

        # Fix grounding
        logger.info('Running grounding map')
        twg = gm.agent_texts_with_grounding(stmts)
        prot_map = gm.protein_map_from_twg(twg)
        gm.default_grounding_map.update(prot_map)
        gmap = gm.GroundingMapper(gm.default_grounding_map)
        stmts = gmap.map_agents(stmts, do_rename=True)

        logger.info('%d Statements after grounding map' % len(stmts))

        # Fix sites
        sm = SiteMapper(default_site_map)
        stmts, _ = sm.map_sites(stmts)

        logger.info('%d Statements with valid sequence' % len(stmts))

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
                    logger.info('Running %s relevance filter' % rel_key)
                    stmts = self._relevance_filter(stmts, [rel_key])
                    logger.info('%s Statements after filter' % len(stmts))

        # Combine duplicates
        logger.info('Preassembling %d Statements' % len(stmts))
        pa = Preassembler(hierarchies, stmts)
        self.unique_stmts = pa.combine_duplicates()
        logger.info('%d unique Statements' % len(self.unique_stmts))

        # Run BeliefEngine on unique statements
        be = BeliefEngine()
        be.set_prior_probs(self.unique_stmts)

        # Build statement hierarchy
        self.toplevel_stmts = pa.combine_related()
        logger.info('%d top-level Statements' % len(self.toplevel_stmts))
        # Run BeliefEngine on hierarchy
        be.set_hierarchy_probs(self.toplevel_stmts)

    def load_prior(self, prior_fname):
        """Load a set of prior statements from a pickle file.

        The prior statements have a special key in the stmts dictionary
        called "prior".

        Parameters
        ----------
        prior_fname : str
            The name of the pickle file containing the prior Statements.
        """
        with open(prior_fname, 'rb') as fh:
            self.stmts['prior'] = pickle.load(fh)

    def get_model_agents(self):
        """Return a list of all Agents from all Statements.

        Returns
        -------
        agents : list[indra.statements.Agent]
           A list of Agents that are in the model.
        """
        model_stmts = self.get_statements_noprior()
        agents = []
        for stmt in model_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.append(a)
        return agents

    def get_prior_agents(self):
        """Return a list of all Agents from the prior Statements.

        Returns
        -------
        agents : list[indra.statements.Agent]
           A list of Agents that are in the prior.
        """
        prior_stmts = self.get_statements_prior()
        agents = []
        for stmt in prior_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.append(a)
        return agents

    def get_statements(self):
        """Return a list of all Statements in a single list.

        Returns
        -------
        stmts : list[indra.statements.Statement]
            A list of all the INDRA Statements in the model.
        """
        stmt_lists = [v for k, v in self.stmts.items()]
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def get_statements_noprior(self):
        """Return a list of all non-prior Statements in a single list.

        Returns
        -------
        stmts : list[indra.statements.Statement]
            A list of all the INDRA Statements in the model (excluding
            the prior).
        """
        stmt_lists = [v for k, v in self.stmts.items() if k != 'prior']
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def get_statements_prior(self):
        """Return a list of all prior Statements in a single list.

        Returns
        -------
        stmts : list[indra.statements.Statement]
            A list of all the INDRA Statements in the prior.
        """
        if self.stmts.get('prior') is not None:
            return self.stmts['prior']
        return []


def _get_agent_comp(agent):
    eh = hierarchies['entity']
    a_ns, a_id = agent.get_grounding()
    if (a_ns is None) or (a_id is None):
        return None
    uri = eh.get_uri(a_ns, a_id)
    comp_id = eh.components.get(uri)
    return comp_id


def _ref_agents_all_filter(stmts_in, ref_agents):
    # If there is no reference, keep everything by default
    if not ref_agents:
        return stmts_in
    stmts_out = []
    # Preprocess reference Agents: make a list of entity hierarchy components
    # that appear in the reference and also a list of reference Agent names
    ref_agent_names = set()
    ref_components = set()
    for a in ref_agents:
        comp_id = _get_agent_comp(a)
        if comp_id is not None:
            ref_components.add(comp_id)
        ref_agent_names.add(a.name)
    # Iterate over every Statement and check if any of its Agents are either
    # in a component appearing in the reference, or match one of the
    # reference Agents that isn't in any of the components.
    for st in stmts_in:
        agents = [a for a in st.agent_list() if a is not None]
        found_all = True
        for st_agent in agents:
            found = False
            comp_id = _get_agent_comp(st_agent)
            if comp_id is None:
                for ref_agent_name in ref_agent_names:
                    if st_agent.name == ref_agent_name:
                        found = True
            elif comp_id in ref_components:
                found = True
            if not found:
                found_all = False
                break
        if found_all:
            stmts_out.append(st)
    return stmts_out


def _ref_agents_one_filter(stmts_in, ref_agents):
    # If there is no reference, keep everything by default
    if not ref_agents:
        return stmts_in
    stmts_out = []
    # Preprocess reference Agents: make a list of entity hierarchy components
    # that appear in the reference and also a list of reference Agent names
    ref_agent_names = set()
    ref_components = set()
    for a in ref_agents:
        comp_id = _get_agent_comp(a)
        if comp_id is not None:
            ref_components.add(comp_id)
        ref_agent_names.add(a.name)

    # Iterate over every Statement and check if any of its Agents are either
    # in a component appearing in the reference, or match one of the
    # reference Agents that isn't in any of the components.
    for st in stmts_in:
        agents = [a for a in st.agent_list() if a is not None]
        found = False
        for st_agent in agents:
            comp_id = _get_agent_comp(st_agent)
            if comp_id is None:
                for ref_agent_name in ref_agent_names:
                    if st_agent.name == ref_agent_name:
                        found = True
                        break
            elif comp_id in ref_components:
                found = True
                break
        if found:
            stmts_out.append(st)
    return stmts_out


def _grounding_filter(stmts_in):
    stmts_out = []
    for st in stmts_in:
        # Check that all agents are grounded
        grounded = True
        for ag in st.agent_list():
            if ag is None:
                continue
            if not ag.db_refs or list(ag.db_refs.keys()) == ['TEXT']:
                grounded = False
                break
        if grounded:
            stmts_out.append(st)
    return stmts_out


def _human_only_filter(stmts_in):
    stmts_out = []
    for st in stmts_in:
        agents = [a for a in st.agent_list() if a is not None]
        non_human = False
        for a in agents:
            hgnc_id = a.db_refs.get('HGNC')
            up_id = a.db_refs.get('UP')
            if not hgnc_id:
                if up_id and not uniprot_client.is_human(up_id):
                    non_human = True
                    break
        if not non_human:
            stmts_out.append(st)
    return stmts_out


def _agent_related(a1, a2):
    if a1.matches(a2) or a1.isa(a2, hierarchies) or a2.isa(a1, hierarchies):
        return True
    return False
