from pysb import *
from pysb.core import SelfExporter
from bel import bel_api
from biopax import biopax_api
from trips import trips_api

SelfExporter.do_export = False

class AgentSet(object):
    """A container for a set of Agents. Wraps a dict of Agent instances."""
    def __init__(self):
        self.agents = {}

    def get_create_agent(self, name):
        """Return agent with given name, creating it if needed."""
        try:
            agent = self.agents[name]
        except KeyError:
            agent = Agent(name)
            self.agents[name] = agent
        return agent

    def iteritems(self):
        return self.agents.iteritems()

    def __getitem__(self, name):
        return self.agents[name]

class Agent(object):
    def __init__(self, name):
        self.name = name
        self.sites = []
        self.site_states = {}
        # The list of site/state configurations that lead to this agent
        # being active (where the agent is currently assumed to have only
        # one type of activity)
        self.activating_mods = []

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
    try:
        default_ic = model.parameters['default_ic']
    except KeyError:
        default_ic = Parameter('default_ic', 100.)
        model.add_component(default_ic)
    # Iterate over all monomers
    for m in model.monomers:
        # Build up monomer pattern dict
        sites_dict = {}
        for site in m.sites:
            if site in m.site_states:
                sites_dict[site] = m.site_states[site][0]
            else:
                sites_dict[site] = None
        mp = m(**sites_dict)
        model.initial(mp, default_ic)


class PysbAssembler(object):
    def __init__(self):
        self.statements = []
        self.agent_set = None

    def add_statements(self, stmts):
        self.statements.extend(stmts)

    def make_model(self, initial_conditions=True):
        model = Model()
        self.agent_set = AgentSet()
        # Collect information about the monomers/self.agent_set from the
        # statements
        for stmt in self.statements:
            stmt.monomers(self.agent_set)
        # Add the monomers to the model based on our AgentSet
        for agent_name, agent in self.agent_set.iteritems():
            m = Monomer(agent_name, agent.sites, agent.site_states)
            model.add_component(m)
        # Iterate over the statements to generate rules
        for stmt in self.statements:
            stmt.assemble(model, self.agent_set)
        # Add initial conditions
        if initial_conditions:
            add_default_initial_conditions(model)
        return model

if __name__ == '__main__':
    pa = PysbAssembler()
    bp = bel_api.process_belrdf('data/RAS_neighborhood.rdf')
    pa.add_statements(bp.belpy_stmts)
    # bp = bel_api.process_ndex_neighborhood("ARAF")
    # pa.add_statements(bp.belpy_stmts)
    # tp = trips_api.process_text("BRAF phosphorylates MEK1 at Ser222")
    # pa.add_statements(tp.belpy_stmts)
    model = pa.make_model()
