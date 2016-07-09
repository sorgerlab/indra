import pickle
import logging
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
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
                self.stmts = pickle.load(open(model_fname, 'rb'))
            except:
                logger.warning('Could not load %s, starting new model.' %
                               model_fname)
                self.stmts = {}
        self.unique_stmts = []
        self.toplevel_stmts = []

    def save(self, model_fname='model.pkl'):
        """Save the state of the IncrementalModel in a pickle file.

        model_fname : Optional[str]
            The name of the pickle file to save the state of the
            IncrementalModel in. Default: model.pkl
        """
        with open(model_fname, 'wb') as fh:
            pickle.dump(self.stmts, fh)

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

        stmts_to_add = range(len(stmts))
        # Filter for grounding
        if 'grounding' in filters:
            for i, stmt in enumerate(stmts):
                agents = [a for a in stmt.agent_list() if a is not None]
                # Check that all agents are grounded
                if any(not a.db_refs for a in agents):
                    stmts_to_add.remove(i)

        if ('prior_all' in filters) or ('prior_one' in filters):
            prior_agents = self.get_prior_agents()
            if prior_agents:
                for i in stmts_to_add:
                    agents = set([a.name for a in stmts[i].agent_list()
                                  if a is not None])
                    if 'prior_all' in filters:
                        if any(not a in prior_agents for a in agents):
                            stmts_to_add.remove(i)
                    if 'prior_one' in filters:
                        if all(not a in prior_agents for a in agents):
                            stmts_to_add.remove(i)

        if ('model_all' in filters) or ('model_one' in filters):
            model_agents = self.get_model_agents()
            if model_agents:
                for i in stmts_to_add:
                    agents = set([a.name for a in stmts[i].agent_list()
                                  if a is not None])
                    if 'model_all' in filters:
                        if any(not a in model_agents for a in agents):
                            stmts_to_add.remove(i)
                    if 'model_one' in filters:
                        if all(not a in model_agents for a in agents):
                            stmts_to_add.remove(i)
        self.stmts[pmid] = [stmts[i] for i in stmts_to_add]

    def preassemble(self):
        """Preassemble the Statements collected in the model.

        Use INDRA's Preassembler on the IncrementalModel and save the
        unique statements and the top level statements in class
        attributes.
        """
        stmts = self.get_statements()
        pa = Preassembler(hierarchies, stmts)
        self.unique_stmts = pa.combine_duplicates()
        self.toplevel_stmts = pa.combine_related()

    def load_prior(self, prior_fname):
        """Load a set of prior statements from a pickle file.

        The prior statements have a special key in the stmts dictionary
        called "prior".
        """
        with open(prior_fname, 'rb') as fh:
            self.stmts['prior'] = pickle.load(fh)

    def get_model_agents(self):
        """Return a list of all Agents from all Statements."""
        model_stmts = self.get_statements()
        agents = set()
        for stmt in model_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add(a.name)
        return agents

    def get_prior_agents(self):
        """Return a list of all Agents from the prior Statements."""
        prior_stmts = self.stmts.get('prior')
        agents = set()
        if prior_stmts is None:
            return agents
        for stmt in prior_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add(a.name)
        return agents

    def get_statements(self):
        """Return a list of all Statements in a single list."""
        stmt_lists = [v for k, v in self.stmts.iteritems()]
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def make_model(self):
        """Assemble a PySB model from all Statements."""
        pa = PysbAssembler()
        pa.add_statements(self.get_statements())
        pa.make_model()
        return pa.model
