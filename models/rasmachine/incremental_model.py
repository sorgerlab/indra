import pickle
from indra.assemblers import PysbAssembler
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import entity_hierarchy as eh
from indra.preassembler.hierarchy_manager import modification_hierarchy as mh

class IncrementalModel(object):
    def __init__(self, fname=None):
        if fname is None:
            self.stmts = {}
        else:
            try:
                self.stmts = pickle.load(open(fname, 'rb'))
            except:
                print 'Could not load %s' % fname
                self.stmts = {}
        self.unique_stmts = []
        self.toplevel_stmts = []

    def save(self, fname='model.pkl'):
        with open(fname, 'wb') as fh:
            pickle.dump(self.stmts, fh)

    def add_statements(self, pmid, stmts, filters=None):
        # Adds statements to the incremental model indexed by
        # PMID. Filters are optional and are a list of filter
        # options. Currently the following options are implemented:
        # - grounding: require that all Agents in statements are grounded
        # - model_one: require that at least one Agent is in the incremental
        #              model
        # - model_all: require that all Agents are in the incremental model
        # - prior_one: require that at least one Agent is in the
        #              prior model
        # - prior_all: require that all Agents are in the prior model
        # Note that model_one -> prior_all are increasingly more restrictive
        # options.

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
        '''Preassemble the statements collected in the model
        and save the unique statements and the top level
        statements'''
        stmts = self.get_statements()
        pa = Preassembler(eh, mh, stmts)
        self.unique_stmts = pa.combine_duplicates()
        self.toplevel_stmts = pa.combine_related()

    def load_prior(self, prior_fname):
        '''Load a set of prior statements from a pickle file'''
        with open(prior_fname, 'rb') as fh:
            self.stmts['prior'] = pickle.load(fh)

    def get_model_agents(self):
        model_stmts = self.get_statements()
        agents = set()
        for stmt in model_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add(a.name)
        return agents

    def get_prior_agents(self):
        prior_stmts = self.stmts.get('prior')
        agents = set()
        if prior_stmts is None:
            return agents
        for stmt in prior_stmts:
            for a in stmt.agent_list():
                if a is not None:
                    agents.add(a.name)
        return agents

    def add_statement(self, pmid, stmt):
        try:
            self.stmts[pmid].append(stmt)
        except KeyError:
            self.stmts[pmid] = [stmt]

    def get_statements(self):
        stmt_lists = [v for k, v in self.stmts.iteritems()]
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def make_model(self):
        pa = PysbAssembler()
        pa.add_statements(self.get_statements())
        pa.make_model()
        return pa.model
