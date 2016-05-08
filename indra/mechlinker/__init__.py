from indra.statements import *

class BaseAgentSet(object):
    """Container for a set of BaseAgents.
    Wraps a dict of BaseAgent instances.
    """
    def __init__(self):
        self.agents = {}

    def get_create_base_agent(self, agent):
        """Return agent with given name, creating it if needed."""
        try:
            base_agent = self.agents[agent.name]
        except KeyError:
            base_agent = BaseAgent(agent.name)
            self.agents[agent.name] = base_agent

        # Handle modification conditions
        for mc in agent.mods:
            base_agent.states.append(mc)

        return base_agent

    def iteritems(self):
        return self.agents.iteritems()

    def __getitem__(self, name):
        return self.agents[name]

class BaseAgent(object):
    def __init__(self, name):
        self.name = name
        self.activities = []
        self.active_states = {}
        self.states = []

    def add_activity(self, activity):
        if activity not in self.activities:
            self.activities.append(activity)
   
    def add_active_state(self, activity, mods):
        if not mods:
            return
        try:
            if not self.hasmod(self.active_states[activity], mods):
                self.active_states[activity].append(mods)
        except KeyError:
            self.active_states[activity] = [mods]

    @staticmethod
    def hasmod(mods_list, mods):
        for ml in mods_list:
            found_ix = []
            for m1 in ml:
                for ix, m2 in enumerate(mods):
                    if m1.equals(m2) and ix not in found_ix:
                        found_ix.append(ix)
                        break
            if len(found_ix) == len(mods):
                return True
        return False

    def __str__(self):
        s = '%s(' % self.name
        if self.activities:
            s += 'activities: %s, ' % self.activities
        for k, v in self.active_states.iteritems():
            s += '%s: %s' % (k, v)
        s += ')'
        return s

    def __repr__(self):
        return self.__str__()

class MechLinker(object):
    def __init__(self, stmts=None):
        if stmts is not None:
            self.statements = stmts
        else:
            self.statements = []
        self.base_agents = BaseAgentSet()

    def add_statements(self, stmts):
        self.statements.extend(stmts)

    def get_activities(self):
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation) or\
                isinstance(stmt, Transphosphorylation) or\
                isinstance(stmt, Autophosphorylation):
                if stmt.enz is not None:
                    enz_base = self.get_base(stmt.enz)
                    enz_base.add_activity('kinase')
                    enz_base.add_active_state('kinase', stmt.enz.mods)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is not None:
                    enz_base = self.get_base(stmt.enz)
                    enz_base.add_activity('phosphatase')
                    enz_base.add_active_state('phosphatase', stmt.enz.mods)
            elif isinstance(stmt, Modification):
                if stmt.enz is not None:
                    enz_base = self.get_base(stmt.enz)
                    enz_base.add_activity('catalytic')
                    enz_base.add_active_state('catalytic', stmt.enz.mods)
            elif isinstance(stmt, SelfModification):
                if stmt.enz is not None:
                    enz_base = self.get_base(stmt.enz)
                    enz_base.add_activity('catalytic')
                    enz_base.add_active_state('catalytic', stmt.enz.mods)
            elif isinstance(stmt, RasGef):
                if stmt.gef is not None:
                    gef_base = self.get_base(stmt.gef)
                    gef_base.add_activity('rasgef')
                    gef_base.add_active_state('rasgef', stmt.gef.mods)
            elif isinstance(stmt, RasGap):
                if stmt.gef is not None:
                    gap_base = self.get_base(stmt.gap)
                    gap_base.add_activity('rasgap')
                    gap_base.add_active_state('rasgap', stmt.gap.mods)
            elif isinstance(stmt, ActivityActivity):
                if stmt.subj is not None:
                    subj_base =\
                        self.get_base(stmt.subj)
                    subj_base.add_activity(stmt.subj_activity)
                    subj_base.add_active_state(stmt.subj_activity,
                                               stmt.subj.mods)
                if stmt.obj is not None:
                    obj_base =\
                        self.get_base(stmt.obj)
                    obj_base.add_activity(stmt.obj_activity)
                    obj_base.add_active_state(stmt.obj_activity,
                                              stmt.obj.mods)
            elif isinstance(stmt, ActiveForm):
                agent_base = self.get_base(stmt.agent)
                agent_base.add_activity(stmt.activity)
                agent_base.add_active_state(stmt.activity, stmt.agent.mods)

    def get_base(self, agent):
        base_agent = self.base_agents.get_create_base_agent(agent)
        return base_agent


if __name__ == '__main__':
    import pickle
    a = pickle.load(open('models/ras_220_genes/ras_genes_results.pkl'))
    ml = MechLinker(a['related2'])
    ml.get_activities() 
