from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import logging
import networkx
import itertools
from copy import deepcopy
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

    def require_active_form(self):
        new_stmts = []
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                if stmt.enz is None:
                    continue
                enz_base = self.get_base(stmt.enz)
                active_forms = enz_base.get_active_forms()
                if not active_forms:
                    new_stmts.append(stmt)
                else:
                    for af in active_forms:
                        new_stmt = deepcopy(stmt)
                        af.apply_to(new_stmt.enz)
                        new_stmts.append(new_stmt)
            elif isinstance(stmt, RegulateAmount) or \
                isinstance(stmt, RegulateActivity):
                if stmt.subj is None:
                    continue
                subj_base = self.get_base(stmt.subj)
                active_forms = subj_base.get_active_forms()
                if not active_forms:
                    new_stmts.append(stmt)
                else:
                    for af in active_forms:
                        new_stmt = deepcopy(stmt)
                        af.apply_to(new_stmt.subj)
                        new_stmts.append(new_stmt)
            else:
                new_stmts.append(stmt)
        self.statements = new_stmts
        return new_stmts


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
                    obj_base = self.get_base(stmt.obj)
                    act_red = \
                        obj_base.get_activity_reduction(stmt.obj_activity)
                    if act_red is not None:
                        stmt.obj_activity = act_red
            elif isinstance(stmt, ActiveForm):
                agent_base = self.get_base(stmt.agent)
                act_red = agent_base.get_activity_reduction(stmt.activity)
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
        linked_stmts = self.infer_active_forms(self.statements)
        linked_stmts += self.infer_modifications(self.statements)
        return linked_stmts

    @staticmethod
    def infer_activations(stmts):
        linked_stmts = []
        act_stmts = get_statement_type(stmts, RegulateActivity)
        af_stmts = get_statement_type(stmts, ActiveForm)
        mod_stmts = get_statement_type(stmts, Modification)
        for af_stmt, mod_stmt in itertools.product(*(af_stmts, mod_stmts)):
            # There has to be an enzyme and the substrate and the
            # agent of the active form have to match
            if mod_stmt.enz is None or \
                (not af_stmt.agent.entity_matches(mod_stmt.sub)):
                continue
            # We now check the modifications to make sure they are consistent
            if not af_stmt.agent.mods:
                continue
            found = False
            for mc in af_stmt.agent.mods:
                if stmt_mod_map.get(mc.mod_type) == mod_stmt.__class__ and \
                    mc.residue == mod_stmt.residue and \
                    mc.position == mod_stmt.position:
                    found = True
            if not found:
                continue
            # Collect evidence
            ev = mod_stmt.evidence
            # Finally, check the polarity of the ActiveForm
            if af_stmt.is_active:
                st = Activation(mod_stmt.enz, mod_stmt.agent, af_stmt.activity,
                                evidence=ev)
            else:
                st = Inhibition(mod_stmt.enz, mod_stmt.agent, af_stmt.activity,
                                evidence=ev)
            linked_stmts.append(LinkedStatement([af_stmt, mod_stmt], st))
        return linked_stmts

    @staticmethod
    def infer_active_forms(stmts):
        # Infer ActiveForm from RegulateActivity + Modification
        linked_stmts = []
        for act_stmt in get_statement_type(stmts, RegulateActivity):
            # TODO: revise the conditions here
            if not (act_stmt.subj.activity is not None and \
                act_stmt.subj.activity.activity_type == 'kinase' and \
                act_stmt.subj.activity.is_active):
                continue
            matching = []
            ev = act_stmt.evidence
            for mod_stmt in get_statement_type(stmts, Modification):
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
        return linked_stmts

    @staticmethod
    def infer_modifications(stmts):
        # Infer indirect Modification from RegulateActivity + ActiveForm
        linked_stmts = []
        for act_stmt in get_statement_type(stmts, RegulateActivity):
            for af_stmt in get_statement_type(stmts, ActiveForm):
                if not af_stmt.agent.entity_matches(act_stmt.obj):
                    continue
                mods = af_stmt.agent.mods
                # Make sure the ActiveForm only involves modified sites
                if af_stmt.agent.mutations or \
                    af_stmt.agent.bound_conditions or \
                    af_stmt.agent.location:
                    continue
                if not af_stmt.agent.mods:
                    continue
                for mod in af_stmt.agent.mods:
                    evs = act_stmt.evidence + af_stmt.evidence
                    for ev in evs:
                        ev.epistemics['direct'] = False
                    if mod.is_modified:
                        mod_type_name = mod.mod_type
                    else:
                        mod_type_name = 'de' + mod.mod_type
                    mod_class = stmt_mod_map.get(mod_type_name)
                    if not mod_class:
                        continue
                    st = mod_class(act_stmt.subj,
                                   act_stmt.obj,
                                   mod.residue, mod.position,
                                   evidence=evs)
                    ls = LinkedStatement([act_stmt, af_stmt], st)
                    linked_stmts.append(ls)
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
    """

    def __init__(self, name):
        self.name = name
        self.activity_types = []
        self.active_states = {}
        self.activity_graph = None
        self.activity_reductions = None

    def get_activity_reduction(self, activity):
        if self.activity_reductions is None:
            self.make_activity_reductions()
        return self.activity_reductions.get(activity)

    def make_activity_reductions(self):
        self.make_activity_graph()
        self.activity_reductions = get_graph_reductions(self.activity_graph)

    def make_activity_graph(self):
        self.activity_graph = networkx.DiGraph()
        for a1, a2 in itertools.combinations(self.activity_types, 2):
            if ah.isa('INDRA', a1, 'INDRA', a2):
                self.activity_graph.add_edge(a2, a1)
            if ah.isa('INDRA', a2, 'INDRA', a1):
                self.activity_graph.add_edge(a1, a2)

    def add_activity(self, activity_type):
        if activity_type not in self.activity_types:
            self.activity_types.append(activity_type)

    def add_active_state(self, activity_type, agent):
        agent_state = AgentState(agent)
        if activity_type in self.active_states:
            self.active_states[activity_type].append(agent_state)
        else:
            self.active_states[activity_type] = [agent_state]

    def get_active_forms(self):
        # TODO: handle inactive states
        # TODO: handle activity types
        if self.active_states:
            states = []
            for k, v in self.active_states.items():
                states += v
            return states
        return None

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

class AgentState(object):
    def __init__(self, agent):
        self.bound_conditions = agent.bound_conditions
        self.mods = agent.mods
        self.mutations = agent.mutations
        self.location = agent.location

    def apply_to(self, agent):
        agent.bound_conditions = self.bound_conditions
        agent.mods = self.mods
        agent.mutations = self.mutations
        agent.location = self.location

    def __repr__(self):
        s = 'AgentState(%s, %s, %s, %s)' % (self.bound_conditions, self.mods,
                                            self.mutations, self.location)
        return s

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

def get_graph_reductions(graph):
    def frontier(g, nd):
        if g.out_degree(nd) == 0:
            return set([nd])
        else:
            frontiers = set()
            for n in g.successors(nd):
                frontiers = frontiers.union(frontier(graph, n))
            return frontiers
    reductions = {}
    nodes_sort = networkx.topological_sort(graph)
    frontiers = [frontier(graph, n) for n in nodes_sort]
    for i, n1 in enumerate(nodes_sort):
        for j, n2 in enumerate(nodes_sort):
            if i > j:
                continue
            if frontiers[i] == frontiers[j]:
                reductions[n1] = n2
    return reductions

stmt_mod_map = {
    'phosphorylation': Phosphorylation,
    'dephosphorylation': Dephosphorylation,
    'autophosphorylation': Autophosphorylation,
    'ubiquitination': Ubiquitination,
    'deubiquitination': Deubiquitination,
    'acetylation': Acetylation,
    'deacetylation': Deacetylation,
    'hydroxylation': Hydroxylation,
    'dehydroxylation': Dehydroxylation,
    'sumoylation': Sumoylation,
    'desumoylation': Desumoylation,
    'glycosylation': Glycosylation,
    'deglycosylation': Deglycosylation,
    'farnesylation': Farnesylation,
    'defarnesylation': Defarnesylation,
    'ribosylation': Ribosylation,
    'deribosylation': Deribosylation,
}
