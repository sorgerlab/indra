__all__ = ['BaseAgentSet', 'BaseAgent']
from pysb import Annotation
from indra.statements import *
from .common import _n
from .sites import states, get_binding_site_name, get_mod_site_name


class BaseAgentSet(object):
    """Container for a dict of BaseAgents with their names as keys."""
    def __init__(self):
        self.agents = {}

    def get_create_base_agent(self, agent):
        """Return base agent with given name, creating it if needed."""
        try:
            base_agent = self.agents[_n(agent.name)]
        except KeyError:
            base_agent = BaseAgent(_n(agent.name))
            self.agents[_n(agent.name)] = base_agent

        # If it's a molecular agent
        if isinstance(agent, Agent):
            # Handle bound conditions
            for bc in agent.bound_conditions:
                bound_base_agent = self.get_create_base_agent(bc.agent)
                bound_base_agent.create_site(get_binding_site_name(agent))
                base_agent.create_site(get_binding_site_name(bc.agent))

            # Handle modification conditions
            for mc in agent.mods:
                base_agent.create_mod_site(mc)

            # Handle mutation conditions
            for mc in agent.mutations:
                res_from = mc.residue_from if mc.residue_from else 'mut'
                res_to = mc.residue_to if mc.residue_to else 'X'
                if mc.position is None:
                    mut_site_name = res_from
                else:
                    mut_site_name = res_from + mc.position

                base_agent.create_site(mut_site_name, states=['WT', res_to])

            # Handle location condition
            if agent.location is not None:
                base_agent.create_site('loc', [_n(agent.location)])

            # Handle activity
            if agent.activity is not None:
                site_name = agent.activity.activity_type
                base_agent.create_site(site_name, ['inactive', 'active'])

        # There might be overwrites here
        for db_name, db_ref in agent.db_refs.items():
            base_agent.db_refs[db_name] = db_ref

        return base_agent

    def items(self):
        """Return items for the set of BaseAgents that this class wraps.
        """
        return self.agents.items()

    def __getitem__(self, name):
        return self.agents[name]


class BaseAgent(object):
    """A BaseAgent aggregates the global properties of an Agent.

    The BaseAgent class aggregates the name, sites, site states, active forms,
    inactive forms and database references of Agents from individual INDRA
    Statements. This allows the PySB Assembler to correctly assemble the
    Monomer signatures in the model.
    """

    def __init__(self, name):
        self.name = name
        self.sites = []
        self.site_states = {}
        self.site_annotations = []
        # The list of site/state configurations that lead to this agent
        # being active (where the agent is currently assumed to have only
        # one type of activity)
        self.active_forms = []
        self.activity_types = []
        self.inactive_forms = []
        self.db_refs = {}

    def create_site(self, site, states=None):
        """Create a new site on an agent if it doesn't already exist."""
        if site not in self.sites:
            self.sites.append(site)
        if states is not None:
            self.site_states.setdefault(site, [])
            try:
                states = list(states)
            except TypeError:
                return
            self.add_site_states(site, states)

    def create_mod_site(self, mc):
        """Create modification site for the BaseAgent from a ModCondition."""
        site_name = get_mod_site_name(mc)
        (unmod_site_state, mod_site_state) = states[mc.mod_type]
        self.create_site(site_name, (unmod_site_state, mod_site_state))
        site_anns = [Annotation((site_name, mod_site_state), mc.mod_type,
                                'is_modification')]
        if mc.residue:
            site_anns.append(Annotation(site_name, mc.residue, 'is_residue'))
        if mc.position:
            site_anns.append(Annotation(site_name, mc.position, 'is_position'))
        self.site_annotations += site_anns

    def add_site_states(self, site, states):
        """Create new states on an agent site if the state doesn't exist."""
        for state in states:
            if state not in self.site_states[site]:
                self.site_states[site].append(state)

    def add_activity_form(self, activity_pattern, is_active):
        """Adds the pattern as an active or inactive form to an Agent.

        Parameters
        ----------
        activity_pattern : dict
            A dictionary of site names and their states.
        is_active : bool
            Is True if the given pattern corresponds to an active state.
        """
        if is_active:
            if activity_pattern not in self.active_forms:
                self.active_forms.append(activity_pattern)
        else:
            if activity_pattern not in self.inactive_forms:
                self.inactive_forms.append(activity_pattern)

    def add_activity_type(self, activity_type):
        """Adds an activity type to an Agent.

        Parameters
        ----------
        activity_type : str
            The type of activity to add such as 'activity', 'kinase',
            'gtpbound'
        """
        if activity_type not in self.activity_types:
            self.activity_types.append(activity_type)

