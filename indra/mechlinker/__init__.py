from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import logging
import itertools
from indra.statements import *
from indra.preassembler.hierarchy_manager import activity_hierarchy as ah

logger = logging.getLogger('mechlinker')


class MechLinker(object):
    """Rewrite the activation pattern of Statements and derive new Statements.

    The mechanism linker (MechLinker) traverses a corpus of Statements and
    uses various inference steps to make the activity types and active
    forms consistent among Statements.
    """
    def __init__(self, stmts=None):
        if stmts is not None:
            self.statements = stmts
        else:
            self.statements = []
        self.base_agents = BaseAgentSet()

    def add_statements(self, stmts):
        """Add statements to the MechLinker.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
        """
        self.statements.extend(stmts)

    def link_statements(self):
        """Run all the steps of mechanism linking and return LinkedStatements.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatement object which contain a tuple with
            a list of source Statements and a Statement that has been derived
            or "linked" from the source Statements.
        """
        self.get_explicit_activities()
        self.reduce_activities()
        linked_stmts = self.link_activations()
        for ls in linked_stmts:
            self.statements.append(ls.inferred_stmt)
        return linked_stmts

    def get_explicit_activities(self):
        """Aggregate all explicit activities given for Agents."""
        for stmt in self.statements:
            agents = stmt.agent_list()
            # Activity types given as ActivityConditions
            for agent in agents:
                if agent is not None and agent.activity is not None:
                    agent_base = self.get_base(agent)
                    agent_base.add_activity(agent.activity.activity_type)
            # Object activities given in RegulateActivity statements
            if isinstance(stmt, RegulateActivity):
                if stmt.obj is not None:
                    obj_base = self.get_base(stmt.obj)
                    obj_base.add_activity(stmt.obj_activity)
            # Activity types given in ActiveForms
            elif isinstance(stmt, ActiveForm):
                agent_base = self.get_base(stmt.agent)
                agent_base.add_activity(stmt.activity)
                agent_base.add_active_state(stmt.activity, stmt.agent)


    def get_implicit_activities(self):
        """Aggregate all implicit activities and active forms of Agents.

        Iterate over self.statements and collect the implied activities
        and active forms of Agents that appear in the Statements. These are
        collected in BaseAgents, one for each type of Agent. For instance,
        there is a single BaseAgent for a type of protein like MAPK1 while
        MAPK1 may appear in a variety of states within Agents used in
        Statements.
        """
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
                    gef_base.add_activity('gef')
                    if stmt.gef.activity is not None:
                        act = stmt.gef.activity.activity_type
                    else:
                        act = 'activity'
                    gef_base.add_active_state(act, stmt.gef.mods)
            elif isinstance(stmt, RasGap):
                if stmt.gap is not None:
                    gap_base = self.get_base(stmt.gap)
                    gap_base.add_activity('gap')
                    if stmt.gap.activity is not None:
                        act = stmt.gap.activity.activity_type
                    else:
                        act = 'activity'
                    gap_base.add_active_state('act', stmt.gap.mods)
            elif isinstance(stmt, RegulateActivity):
                if stmt.subj is not None:
                    subj_base = self.get_base(stmt.subj)
                    subj_base.add_activity(stmt.j)
            elif isinstance(stmt, ActiveForm):
                agent_base = self.get_base(stmt.agent)
                agent_base.add_activity(stmt.activity)
                agent_base.add_active_state(stmt.activity, stmt.agent.mods)

    def get_base(self, agent):
        """Return the BaseAgent corresponding to an Agent.

        Parameters
        ----------
        agent : indra.statements.Agent

        Returns
        -------
        base_agent : indra.mechlinker.BaseAgent
        """
        base_agent = self.base_agents.get_create_base_agent(agent)
        return base_agent

    def reduce_activities(self):
        """Rewrite the activity types referenced in Statements for consistency.

        Activity types are reduced to the most specific form whenever possible.
        For instance, if 'kinase' is the only specific activity type known
        for the BaseAgent of BRAF, its generic 'activity' forms are rewritten
        to 'kinase'.
        """
        for stmt in self.statements:
            agents = stmt.agent_list()
            for agent in agents:
                if agent is not None and agent.activity is not None:
                    agent_base = self.get_base(agent)
                    act_red = agent_base.get_activity_reduction(
                                                agent.activity.activity_type)
                    if act_red is not None:
                        agent.activity.activity_type = act_red
            if isinstance(stmt, RegulateActivity):
                if stmt.obj is not None:
                    obj_base =\
                        self.get_base(stmt.obj)
                    act_red = \
                        obj_base.get_activity_reduction(stmt.obj_activity)
                    if act_red is not None:
                        stmt.obj_activity = act_red
            elif isinstance(stmt, ActiveForm):
                agent_base = self.get_base(stmt.agent)
                act_red = \
                    agent_base.get_activity_reduction(stmt.activity)
                if act_red is not None:
                    stmt.activity = act_red

    def link_activations(self):
        """Link Activation/Inhibition, Modification and ActiveForm Statements.

        Finds equivalences between Activation/Inhibition, Modification and
        ActiveForm Statements and derives LinkedStatements from a list of
        source statements.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatement object which contain a tuple with
            a list of source Statements and a Statement that has been derived
            or "linked" from the source Statements.
        """
        linked_stmts = []
        for act_stmt in get_statement_type(self.statements, RegulateActivity):
            # Infer ActiveForm from RegulateActivity + Modification
            if act_stmt.subj.activity is not None and \
                act_stmt.subj.activity.activity_type == 'kinase' and \
                act_stmt.subj.activity.is_active:
                matching = []
                ev = act_stmt.evidence
                for mod_stmt in get_statement_type(self.statements,
                                                   Modification):
                    if mod_stmt.enz is not None:
                        if mod_stmt.enz.entity_matches(act_stmt.subj) and \
                            mod_stmt.sub.entity_matches(act_stmt.obj):
                            matching.append(mod_stmt)
                            ev.extend(mod_stmt.evidence)
                if not matching:
                    continue
                mods = []
                for mod_stmt in matching:
                    mod_type_name = mod_stmt.__class__.__name__.lower()
                    if isinstance(mod_stmt, AddModification):
                        is_modified = True
                    else:
                        is_modified = False
                        mod_type_name = mod_type_name[2:]
                    mc = ModCondition(mod_type_name, mod_stmt.residue,
                                      mod_stmt.position, is_modified)
                    mods.append(mc)
                source_stmts = [act_stmt] + [m for m in matching]
                st = ActiveForm(Agent(act_stmt.obj.name, mods=mods,
                                      db_refs=act_stmt.obj.db_refs),
                                act_stmt.obj_activity, act_stmt.is_activation,
                                evidence=ev)
                linked_stmts.append(LinkedStatement(source_stmts, st))
                logger.info('inferred: %s' % st)
        # Infer indirect Phosphorylation from ActAct + ActiveForm
        for act_stmt in get_statement_type(self.statements, RegulateActivity):
            for af_stmt in get_statement_type(self.statements, ActiveForm):
                if not af_stmt.agent.name == act_stmt.obj.name:
                    continue
                mods = af_stmt.agent.mods
                # Make sure the ActiveForm only involves phosphorylated
                # sites
                if not af_stmt.agent.mutations and \
                    not af_stmt.agent.bound_conditions and \
                    all([m.mod_type == 'phosphorylation' for m in mods]) and \
                    all([m.is_modified for m in mods]):
                    for m in mods:
                        ev = act_stmt.evidence
                        ev[0].epistemics['direct'] = False
                        if (act_stmt.is_activation and af_stmt.is_active) or \
                            (not act_stmt.is_activation and not
                             af_stmt.is_active):
                            st = Phosphorylation(act_stmt.subj,
                                                 act_stmt.obj,
                                                 m.residue, m.position,
                                                 evidence=ev)
                        elif (not act_stmt.is_activation and af_stmt.is_active) or \
                              (act_stmt.is_activation and not af_stmt.is_active):
                            st = Dephosphorylation(act_stmt.subj,
                                                   act_stmt.obj,
                                                   m.residue, m.position,
                                                   evidence=ev)
                        linked_stmts.append(LinkedStatement([act_stmt, af_stmt], st))
                        logger.info('inferred: %s' % st)
        return linked_stmts

class BaseAgentSet(object):
    """Container for a set of BaseAgents.

    This class wraps a dict of BaseAgent instance and can be used to get and
    set BaseAgents.
    """
    def __init__(self):
        self.agents = {}

    def get_create_base_agent(self, agent):
        """Return BaseAgent from an Agent, creating it if needed.

        Parameters
        ----------
        agent : indra.statements.Agent

        Returns
        -------
        base_agent : indra.mechlinker.BaseAgent
        """
        try:
            base_agent = self.agents[agent.name]
        except KeyError:
            base_agent = BaseAgent(agent.name)
            self.agents[agent.name] = base_agent

        # Handle modification conditions
        for mc in agent.mods:
            base_agent.states.append(mc)

        return base_agent

    def keys(self):
        return self.agents.keys()

    def items(self):
        return self.agents.items()

    def __getitem__(self, name):
        return self.agents[name]

@python_2_unicode_compatible
class BaseAgent(object):
    """Represents all activity types and active forms of an Agent.

    Parameters
    ----------
    name : str
        The name of the BaseAgent
    activities : list[str]
        A list of activity types that the Agent has
    active_states : dict
        A dict of activity types and their associated Agent states
    activity_reductions : dict
        A dict of activity types and the type they are reduced to by inference.
    states : list[indra.statements.ModCondition]
        A list of ModConditions that the associated Agent can have
    """
    def AgentState(object):
        def __init__(agent):
            self.bound_conditions = agent.bound_conditions
            self.mods = agent.mods
            self.mutations = agent.mutations
            self.location = agent.location

    def __init__(self, name):
        self.name = name
        self.activity_types = []
        self.active_states = {}
        self.states = []
        self.activity_graph = None
        self.activity_reductions = None

    def get_activity_reduction(self, activity):
        self.make_activity_reductions()
        return self.activity_reductions.get(activity)

    def make_activity_reductions(self):
        self.make_activity_graph()
        self.activity_reductions = get_graph_reductions(self.activity_graph)

    def make_activity_graph(self):
        self.activity_graph  = []
        for a1, a2 in itertools.combinations(self.activity_types, 2):
            if ah.isa('INDRA', a1, 'INDRA', a2):
                self.activity_graph.append((a1, a2))
            if ah.isa('INDRA', a2, 'INDRA', a1):
                self.activity_graph.append((a2, a1))

    def add_activity(self, activity_type):
        if activity_type not in self.activity_types:
            self.activity_types.append(activity_type)

    def add_active_state(self, activity_type, agent):
        agent_state = AgentState(agent)
        self.active_states[activity_type].append(agent_state)

    def __str__(self):
        s = '%s(' % self.name
        if self.activity_types:
            s += 'activity_types: %s, ' % self.activity_types
        for k, v in self.active_states.items():
            s += '%s: %s' % (k, v)
        s += ')'
        return s

    def __repr__(self):
        return str(self)

@python_2_unicode_compatible
class LinkedStatement(object):
    """A tuple containing a list of source Statements and an inferred Statement.

    The list of source Statements are the basis for the inferred Statement.

    Parameters
    ----------
    source_stmts : list[indra.statements.Statement]
        A list of source Statements
    inferred_stmts : indra.statements.Statement
        A Statement that was inferred from the source Statements.
    """
    def __init__(self, source_stmts, inferred_stmt):
        self.source_stmts = source_stmts
        self.inferred_stmt = inferred_stmt

    def __str__(self):
        source_str = ', '.join([str(st) for st in self.source_stmts])
        inferred_str = str(self.inferred_stmt)
        s = 'LinkedStatement((%s), %s)' % (source_str, inferred_str)
        return s

    def __repr__(self):
        return str(self)

def get_statement_type(stmts, stmt_type):
    return [st for st in stmts if isinstance(st, stmt_type)]

def get_graph_reductions(edges):
    reductions = {}
    nodes = set()
    for s, t in edges:
        nodes.add(s)
        nodes.add(t)
    reverse_edges = {n:[] for n in nodes}
    for s, t in edges:
        reverse_edges[t].append(s)
    for n in nodes:
        next_nodes = reverse_edges[n]
        reduced_to = None
        while len(next_nodes) == 1:
            reduced_to = next_nodes[0]
            next_nodes = reverse_edges[next_nodes[0]]
        if reduced_to is not None:
            reductions[n] = reduced_to
    return reductions
