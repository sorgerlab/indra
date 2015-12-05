from pysb import Model, Monomer, Parameter, Annotation
from pysb.core import SelfExporter
import pysb.export
from bel import bel_api
from biopax import biopax_api
from trips import trips_api
from indra.statements import Agent

SelfExporter.do_export = False

class BaseAgentSet(object):
    """A container for a set of BaseAgents. Wraps a dict of BaseAgent instances."""
    def __init__(self):
        self.agents = {}

    def get_create_base_agent(self, agent):
        """Return agent with given name, creating it if needed."""
        try:
            base_agent = self.agents[agent.name]
        except KeyError:
            base_agent = BaseAgent(agent.name)
            self.agents[agent.name] = base_agent
        
        if agent.bound_to:
            bound_to = self.get_create_base_agent(Agent(agent.bound_to))
            bound_to.create_site(agent.name)
            base_agent.create_site(agent.bound_to)
       
        # There might be overwrites here
        for db_name, db_ref in agent.db_refs.iteritems():
            base_agent.db_refs[db_name] = db_ref

        return base_agent

    def iteritems(self):
        return self.agents.iteritems()

    def __getitem__(self, name):
        return self.agents[name]

class BaseAgent(object):
    def __init__(self, name):
        self.name = name
        self.sites = []
        self.site_states = {}
        # The list of site/state configurations that lead to this agent
        # being active (where the agent is currently assumed to have only
        # one type of activity)
        self.activating_mods = []
        self.db_refs = {}

    def create_site(self, site, states=None):
        """Create a new site on an agent if it doesn't already exist"""
        if site not in self.sites:
            self.sites.append(site)
        if states is not None:
            self.site_states.setdefault(site, [])
            try:
                states = list(states)
            except TypeError:
                return
            self.add_site_states(site, states)

    def add_site_states(self, site, states):
        """Create new states on a agent site if the site doesn't exist"""
        for state in states:
            if state not in self.site_states[site]:
                self.site_states[site].append(state)

    def add_activating_modification(self, activity_pattern):
        self.activating_mods.append(activity_pattern)

def add_default_initial_conditions(model):
    # Iterate over all monomers
    for m in model.monomers:
        set_base_initial_condition(model, m, 100.0)

def set_base_initial_condition(model, monomer, value):
    # Build up monomer pattern dict
    sites_dict = {}
    for site in monomer.sites:
        if site in monomer.site_states:
            sites_dict[site] = monomer.site_states[site][0]
        else:
            sites_dict[site] = None
    mp = monomer(**sites_dict)
    pname = monomer.name + '_0'
    try:
        p = model.parameters[pname]
        p.value = value
    except KeyError:
        p = Parameter(pname, value)
        model.add_component(p)
        model.initial(mp, p)
    
def get_annotation(component, db_name, db_ref):
    '''
    Construct Annotation following format guidelines 
    given at http://identifiers.org/.
    '''
    url = 'http://identifiers.org/'
    subj = component
    if db_name == 'UP':
        obj = url + 'uniprot/%s' % db_ref
        pred = 'is'
    elif db_name == 'HGNC':
        obj = url + 'hgnc/HGNC:%s' % db_ref
        pred = 'is'
    elif db_name == 'XFAM' and db_ref.startswith('PF'):
        obj = url + 'pfam/%s' % db_ref
        pred = 'is'
    elif db_name == 'CHEBI':
        obj = url + 'chebi/CHEBI:%s' % db_ref
        pred = 'is'
    else:
        return None
    return Annotation(subj, obj, pred)

class PysbAssembler(object):
    def __init__(self):
        self.statements = []
        self.agent_set = None
        self.model = None

    def statement_exists(self, stmt):
        for s in self.statements:
            if stmt == s:
                return True
        return False

    def add_statements(self, stmts):
        for stmt in stmts:
            if not self.statement_exists(stmt):
                self.statements.append(stmt)

    def make_model(self, initial_conditions=True, policies=None):
        model = Model()
        # Keep track of which policies we're using
        self.policies = policies
        self.agent_set = BaseAgentSet()
        # Collect information about the monomers/self.agent_set from the
        # statements
        for stmt in self.statements:
            stmt.monomers(self.agent_set, policies=policies)
        # Add the monomers to the model based on our BaseAgentSet
        for agent_name, agent in self.agent_set.iteritems():
            m = Monomer(agent_name, agent.sites, agent.site_states)
            model.add_component(m)
            for db_name, db_ref in agent.db_refs.iteritems():
                a = get_annotation(m, db_name, db_ref)
                if a is not None:
                    model.add_annotation(a)
        # Iterate over the statements to generate rules
        for stmt in self.statements:
            stmt.assemble(model, self.agent_set, policies=policies)
        # Add initial conditions
        if initial_conditions:
            add_default_initial_conditions(model)
        self.model = model
        return self.model
    
    def print_model(self, fname='pysb_model.py'):
        if self.model is not None:
            with open(fname, 'wt') as fh:
                fh.write(pysb.export.export(self.model, 'pysb_flat'))

if __name__ == '__main__':
    pa = PysbAssembler()
    bp = bel_api.process_belrdf('data/RAS_neighborhood.rdf')
    pa.add_statements(bp.statements)
    # bp = bel_api.process_ndex_neighborhood("ARAF")
    # pa.add_statements(bp.statements)
    # tp = trips_api.process_text("BRAF phosphorylates MEK1 at Ser222")
    # pa.add_statements(tp.statements)
    model = pa.make_model()
