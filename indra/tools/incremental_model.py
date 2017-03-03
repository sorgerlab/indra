from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import logging
from indra.statements import Agent
import indra.tools.assemble_corpus as ac
from indra.databases import hgnc_client
from indra.preassembler.hierarchy_manager import hierarchies

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
    assembled_stmts : list[indra.statements.Statement]
        A list of INDRA Statements after assembly.
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
        self.prior_genes = []
        self.assembled_stmts = []

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

    def add_statements(self, pmid, stmts):
        """Add INDRA Statements to the incremental model indexed by PMID.

        Parameters
        ----------
        pmid : str
            The PMID of the paper from which statements were extracted.
        stmts : list[indra.statements.Statement]
            A list of INDRA Statements to be added to the model.
        """
        if pmid not in self.stmts:
            self.stmts[pmid] = stmts
        else:
            self.stmts[pmid] += stmts

    def _relevance_filter(self, stmts, filters=None):
        if filters is None:
            return stmts
        logger.info('Running relevance filter on %d statements' % len(stmts))
        prior_agents = get_gene_agents(self.prior_genes)
        if 'prior_all' in filters:
            stmts = _ref_agents_all_filter(stmts, prior_agents)
        elif 'prior_one' in filters:
            stmts = _ref_agents_one_filter(stmts, prior_agents)
        logger.info('%d statements after relevance filter' % len(stmts))
        return stmts

    def preassemble(self, filters=None):
        """Preassemble the Statements collected in the model.

        Use INDRA's GroundingMapper, Preassembler and BeliefEngine
        on the IncrementalModel and save the unique statements and
        the top level statements in class attributes.

        Currently the following filter options are implemented:
        - grounding: require that all Agents in statements are grounded
        - human_only: require that all proteins are human proteins
        - prior_one: require that at least one Agent is in the prior model
        - prior_all: require that all Agents are in the prior model

        Parameters
        ----------
        filters : Optional[list[str]]
            A list of filter options to apply when choosing the statements.
            See description above for more details. Default: None
        """
        stmts = self.get_statements()

        # Fix grounding
        stmts = ac.map_grounding(stmts)

        if filters and ('grounding' in filters):
            stmts = ac.filter_grounded_only(stmts)

        # Fix sites
        stmts = ac.map_sequence(stmts)

        if filters and 'human_only' in filters:
            stmts = ac.filter_human_only(stmts)

        # Run preassembly
        stmts = ac.run_preassembly(stmts, return_toplevel=False)

        # Run relevance filter
        stmts = self._relevance_filter(stmts, filters)

        # Save Statements
        self.assembled_stmts = stmts

    def load_prior(self, prior_fname):
        """Load a set of prior statements from a pickle file.

        The prior statements have a special key in the stmts dictionary
        called "prior".

        Parameters
        ----------
        prior_fname : str
            The name of the pickle file containing the prior Statements.
        """
        self.stmts['prior'] = ac.load_statements(prior_fname)

    def get_model_agents(self):
        """Return a list of all Agents from all Statements.

        Returns
        -------
        agents : list[indra.statements.Agent]
           A list of Agents that are in the model.
        """
        model_stmts = self.get_statements()
        agents = []
        for stmt in model_stmts:
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

def get_gene_agents(gene_names):
    agents = []
    for gn in gene_names:
        hgnc_id = hgnc_client.get_hgnc_id(gn)
        if not hgnc_id:
            logger.warning('Invalid HGNC gene symbol: %s' % gn)
            continue
        db_refs = {'HGNC': hgnc_id}
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        if up_id:
            db_refs['UP'] = up_id
        agent = Agent(gn, db_refs=db_refs)
        agents.append(agent)
    return agents

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

def _agent_related(a1, a2):
    if a1.matches(a2) or a1.isa(a2, hierarchies) or a2.isa(a1, hierarchies):
        return True
    return False
