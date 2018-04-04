from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import uuid
import logging
import networkx
import itertools
from indra.util import fast_deepcopy
from indra.statements import *
from indra.preassembler.hierarchy_manager import hierarchies

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
            A list of Statements to add.
        """
        self.statements.extend(stmts)

    def gather_explicit_activities(self):
        """Aggregate all explicit activities and active forms of Agents.

        This function iterates over self.statements and extracts explicitly
        stated activity types and active forms for Agents.
        """
        for stmt in self.statements:
            agents = stmt.agent_list()
            # Activity types given as ActivityConditions
            for agent in agents:
                if agent is not None and agent.activity is not None:
                    agent_base = self._get_base(agent)
                    agent_base.add_activity(agent.activity.activity_type)
            # Object activities given in RegulateActivity statements
            if isinstance(stmt, RegulateActivity):
                if stmt.obj is not None:
                    obj_base = self._get_base(stmt.obj)
                    obj_base.add_activity(stmt.obj_activity)
            # Activity types given in ActiveForms
            elif isinstance(stmt, ActiveForm):
                agent_base = self._get_base(stmt.agent)
                agent_base.add_activity(stmt.activity)
                if stmt.is_active:
                    agent_base.add_active_state(stmt.activity, stmt.agent,
                                                stmt.evidence)
                else:
                    agent_base.add_inactive_state(stmt.activity, stmt.agent,
                                                  stmt.evidence)

    def gather_implicit_activities(self):
        """Aggregate all implicit activities and active forms of Agents.

        Iterate over self.statements and collect the implied activities
        and active forms of Agents that appear in the Statements.

        Note that using this function to collect implied Agent activities can
        be risky. Assume, for instance, that a Statement from a reading
        system states that EGF bound to EGFR phosphorylates ERK. This would
        be interpreted as implicit evidence for the EGFR-bound form of EGF
        to have 'kinase' activity, which is clearly incorrect.

        In contrast the alternative pair of this function:
        gather_explicit_activities collects only explicitly stated activities.
        """
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation) or \
                isinstance(stmt, Transphosphorylation) or \
                isinstance(stmt, Autophosphorylation):
                if stmt.enz is not None:
                    enz_base = self._get_base(stmt.enz)
                    enz_base.add_activity('kinase')
                    enz_base.add_active_state('kinase', stmt.enz.mods)
            elif isinstance(stmt, Dephosphorylation):
                if stmt.enz is not None:
                    enz_base = self._get_base(stmt.enz)
                    enz_base.add_activity('phosphatase')
                    enz_base.add_active_state('phosphatase', stmt.enz.mods)
            elif isinstance(stmt, Modification):
                if stmt.enz is not None:
                    enz_base = self._get_base(stmt.enz)
                    enz_base.add_activity('catalytic')
                    enz_base.add_active_state('catalytic', stmt.enz.mods)
            elif isinstance(stmt, SelfModification):
                if stmt.enz is not None:
                    enz_base = self._get_base(stmt.enz)
                    enz_base.add_activity('catalytic')
                    enz_base.add_active_state('catalytic', stmt.enz.mods)
            elif isinstance(stmt, Gef):
                if stmt.gef is not None:
                    gef_base = self._get_base(stmt.gef)
                    gef_base.add_activity('gef')
                    if stmt.gef.activity is not None:
                        act = stmt.gef.activity.activity_type
                    else:
                        act = 'activity'
                    gef_base.add_active_state(act, stmt.gef.mods)
            elif isinstance(stmt, Gap):
                if stmt.gap is not None:
                    gap_base = self._get_base(stmt.gap)
                    gap_base.add_activity('gap')
                    if stmt.gap.activity is not None:
                        act = stmt.gap.activity.activity_type
                    else:
                        act = 'activity'
                    gap_base.add_active_state('act', stmt.gap.mods)
            elif isinstance(stmt, RegulateActivity):
                if stmt.subj is not None:
                    subj_base = self._get_base(stmt.subj)
                    subj_base.add_activity(stmt.j)

    def gather_modifications(self):
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                sub_base = self._get_base(stmt.sub)
                pol = isinstance(stmt, AddModification)
                mod_type = modclass_to_modtype[stmt.__class__]
                if not pol:
                    mod_type = modtype_to_inverse[mod_type]
                mc = ModCondition(mod_type, stmt.residue, stmt.position, pol)
                sub_base.add_modification(mc)
            for agent in stmt.agent_list():
                if agent is not None:
                    agent_base = self._get_base(agent)
                    for mc in agent.mods:
                        agent_base.add_modification(mc)

    def reduce_modifications(self):
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                pol = isinstance(stmt, AddModification)
                mod_type = modclass_to_modtype[stmt.__class__]
                if not pol:
                    mod_type = modtype_to_inverse[mod_type]
                mc = ModCondition(mod_type, stmt.residue, stmt.position, pol)
                sub_base = self._get_base(stmt.sub)
                mc_red = sub_base.get_modification_reduction(mc)
                stmt.residue = mc_red.residue
                stmt.position = mc_red.position
            agents = stmt.agent_list()
            for agent in agents:
                if agent is not None and agent.mods:
                    agent_base = self._get_base(agent)
                    for i, mc in enumerate(agent.mods):
                        mc_red = agent_base.get_modification_reduction(mc)
                        agent.mods[i] = mc_red


    def require_active_forms(self):
        """Rewrites Statements with Agents' active forms in active positions.

        As an example, the enzyme in a Modification Statement can be expected
        to be in an active state. Similarly, subjects of RegulateAmount and
        RegulateActivity Statements can be expected to be in an active form.
        This function takes the collected active states of Agents in their
        corresponding BaseAgents and then rewrites other Statements to apply
        the active Agent states to them.

        Returns
        -------
        new_stmts : list[indra.statements.Statement]
            A list of Statements which includes the newly rewritten Statements.
            This list is also set as the internal Statement list of the
            MechLinker.
        """
        logger.info('Setting required active forms on %d statements...' %
                    len(self.statements))
        new_stmts = []
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                if stmt.enz is None:
                    new_stmts.append(stmt)
                    continue
                enz_base = self._get_base(stmt.enz)
                active_forms = enz_base.get_active_forms()
                if not active_forms:
                    new_stmts.append(stmt)
                else:
                    for af in active_forms:
                        new_stmt = fast_deepcopy(stmt)
                        new_stmt.uuid = str(uuid.uuid4())
                        evs = af.apply_to(new_stmt.enz)
                        new_stmt.partial_evidence = evs
                        new_stmts.append(new_stmt)
            elif isinstance(stmt, RegulateAmount) or \
                isinstance(stmt, RegulateActivity):
                if stmt.subj is None:
                    new_stmts.append(stmt)
                    continue
                subj_base = self._get_base(stmt.subj)
                active_forms = subj_base.get_active_forms()
                if not active_forms:
                    new_stmts.append(stmt)
                else:
                    for af in active_forms:
                        new_stmt = fast_deepcopy(stmt)
                        new_stmt.uuid = str(uuid.uuid4())
                        evs = af.apply_to(new_stmt.subj)
                        new_stmt.partial_evidence = evs
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
                    agent_base = self._get_base(agent)
                    act_red = agent_base.get_activity_reduction(
                                                agent.activity.activity_type)
                    if act_red is not None:
                        agent.activity.activity_type = act_red
            if isinstance(stmt, RegulateActivity):
                if stmt.obj is not None:
                    obj_base = self._get_base(stmt.obj)
                    act_red = \
                        obj_base.get_activity_reduction(stmt.obj_activity)
                    if act_red is not None:
                        stmt.obj_activity = act_red
            elif isinstance(stmt, ActiveForm):
                agent_base = self._get_base(stmt.agent)
                act_red = agent_base.get_activity_reduction(stmt.activity)
                if act_red is not None:
                    stmt.activity = act_red

    @staticmethod
    def infer_complexes(stmts):
        """Return inferred Complex from Statements implying physical interaction.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of Statements to infer Complexes from.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatements representing the inferred Statements.
        """
        interact_stmts = _get_statements_by_type(stmts, Modification)
        linked_stmts = []
        for mstmt in interact_stmts:
            if mstmt.enz is None:
                continue
            st = Complex([mstmt.enz, mstmt.sub], evidence=mstmt.evidence)
            linked_stmts.append(st)
        return linked_stmts

    @staticmethod
    def infer_activations(stmts):
        """Return inferred RegulateActivity from Modification + ActiveForm.

        This function looks for combinations of Modification and ActiveForm
        Statements and infers Activation/Inhibition Statements from them.
        For example, if we know that A phosphorylates B, and the
        phosphorylated form of B is active, then we can infer that
        A activates B. This can also be viewed as having "explained" a given
        Activation/Inhibition Statement with a combination of more mechanistic
        Modification + ActiveForm Statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of Statements to infer RegulateActivity from.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatements representing the inferred Statements.
        """
        linked_stmts = []
        af_stmts = _get_statements_by_type(stmts, ActiveForm)
        mod_stmts = _get_statements_by_type(stmts, Modification)
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
                if mc.mod_type == modclass_to_modtype[mod_stmt.__class__] and \
                    mc.residue == mod_stmt.residue and \
                    mc.position == mod_stmt.position:
                    found = True
            if not found:
                continue
            # Collect evidence
            ev = mod_stmt.evidence
            # Finally, check the polarity of the ActiveForm
            if af_stmt.is_active:
                st = Activation(mod_stmt.enz, mod_stmt.sub, af_stmt.activity,
                                evidence=ev)
            else:
                st = Inhibition(mod_stmt.enz, mod_stmt.sub, af_stmt.activity,
                                evidence=ev)
            linked_stmts.append(LinkedStatement([af_stmt, mod_stmt], st))
        return linked_stmts

    @staticmethod
    def infer_active_forms(stmts):
        """Return inferred ActiveForm from RegulateActivity + Modification.

        This function looks for combinations of Activation/Inhibition
        Statements and Modification Statements, and infers an ActiveForm
        from them. For example, if we know that A activates B and
        A phosphorylates B, then we can infer that the phosphorylated form
        of B is active.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of Statements to infer ActiveForms from.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatements representing the inferred Statements.
        """
        linked_stmts = []
        for act_stmt in _get_statements_by_type(stmts, RegulateActivity):
            # TODO: revise the conditions here
            if not (act_stmt.subj.activity is not None and
                act_stmt.subj.activity.activity_type == 'kinase' and
                act_stmt.subj.activity.is_active):
                continue
            matching = []
            ev = act_stmt.evidence
            for mod_stmt in _get_statements_by_type(stmts, Modification):
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
        """Return inferred Modification from RegulateActivity + ActiveForm.

        This function looks for combinations of Activation/Inhibition Statements
        and ActiveForm Statements that imply a Modification Statement.
        For example, if we know that A activates B, and phosphorylated B is
        active, then we can infer that A leads to the phosphorylation of B.
        An additional requirement when making this assumption is that the
        activity of B should only be dependent on the modified state and not
        other context - otherwise the inferred Modification is not necessarily
        warranted.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of Statements to infer Modifications from.

        Returns
        -------
        linked_stmts : list[indra.mechlinker.LinkedStatement]
            A list of LinkedStatements representing the inferred Statements.
        """
        linked_stmts = []
        for act_stmt in _get_statements_by_type(stmts, RegulateActivity):
            for af_stmt in _get_statements_by_type(stmts, ActiveForm):
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
                        mod_type_name = modtype_to_inverse[mod.mod_type]
                    mod_class = modtype_to_modclass[mod_type_name]
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

    def replace_complexes(self, linked_stmts=None):
        """Remove Complex Statements that can be inferred out.

        This function iterates over self.statements and looks for Complex
        Statements that either match or are refined by inferred Complex
        Statements that were linked (provided as the linked_stmts argument).
        It removes Complex Statements from self.statements that can be
        explained by the linked statements.

        Parameters
        ----------
        linked_stmts : Optional[list[indra.mechlinker.LinkedStatement]]
            A list of linked statements, optionally passed from outside.
            If None is passed, the MechLinker runs self.infer_complexes to
            infer Complexes and obtain a list of LinkedStatements that are
            then used for removing existing Complexes in self.statements.
        """
        if linked_stmts is None:
            linked_stmts = self.infer_complexes(self.statements)
        new_stmts = []
        for stmt in self.statements:
            if not isinstance(stmt, Complex):
                new_stmts.append(stmt)
                continue
            found = False
            for linked_stmt in linked_stmts:
                if linked_stmt.refinement_of(stmt, hierarchies):
                    found = True
            if not found:
                new_stmts.append(stmt)
            else:
                logger.info('Removing complex: %s' % stmt)
        self.statements = new_stmts

    def replace_activations(self, linked_stmts=None):
        """Remove RegulateActivity Statements that can be inferred out.

        This function iterates over self.statements and looks for
        RegulateActivity Statements that either match or are refined by
        inferred RegulateActivity Statements that were linked
        (provided as the linked_stmts argument).
        It removes RegulateActivity Statements from self.statements that can be
        explained by the linked statements.

        Parameters
        ----------
        linked_stmts : Optional[list[indra.mechlinker.LinkedStatement]]
            A list of linked statements, optionally passed from outside.
            If None is passed, the MechLinker runs self.infer_activations to
            infer RegulateActivities and obtain a list of LinkedStatements
            that are then used for removing existing Complexes
            in self.statements.
        """
        if linked_stmts is None:
            linked_stmts = self.infer_activations(self.statements)
        new_stmts = []
        for stmt in self.statements:
            if not isinstance(stmt, RegulateActivity):
                new_stmts.append(stmt)
                continue
            found = False
            for linked_stmt in linked_stmts:
                inferred_stmt = linked_stmt.inferred_stmt
                if stmt.is_activation == inferred_stmt.is_activation and \
                    stmt.subj.entity_matches(inferred_stmt.subj) and \
                    stmt.obj.entity_matches(inferred_stmt.obj):
                        found = True
            if not found:
                new_stmts.append(stmt)
            else:
                logger.info('Removing regulate activity: %s' % stmt)
        self.statements = new_stmts

    def _get_base(self, agent):
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
    activity_types : list[str]
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
        self.inactive_states = {}
        self.activity_graph = None
        self.activity_reductions = None
        self.modification_reductions = None
        self.modifications = []

    def get_activity_reduction(self, activity):
        if self.activity_reductions is None:
            self._make_activity_reductions()
        return self.activity_reductions.get(activity)

    def _make_activity_reductions(self):
        self._make_activity_graph()
        self.activity_reductions = _get_graph_reductions(self.activity_graph)

    def _make_activity_graph(self):
        self.activity_graph = networkx.DiGraph()
        for a1, a2 in itertools.combinations(self.activity_types, 2):
            if hierarchies['activity'].isa('INDRA_ACTIVITIES', a1,
                                           'INDRA_ACTIVITIES', a2):
                self.activity_graph.add_edge(a2, a1)
            if hierarchies['activity'].isa('INDRA_ACTIVITIES', a2,
                                           'INDRA_ACTIVITIES', a1):
                self.activity_graph.add_edge(a1, a2)

    def get_modification_reduction(self, mc):
        if self.modification_reductions is None:
            self._make_modification_reductions()
        mc_red_tuple = self.modification_reductions.get(_mc_tuple(mc))
        # This handles the case where there was no reduction
        if not mc_red_tuple:
            return mc
        mc = ModCondition(*(list(mc_red_tuple) + [mc.is_modified]))
        return mc

    def _make_modification_reductions(self):
        self._make_modification_graph()
        self.modification_reductions = \
            _get_graph_reductions(self.modification_graph)

    def _make_modification_graph(self):
        self.modification_graph = networkx.DiGraph()
        for m1, m2 in itertools.combinations(self.modifications, 2):
            if m1.refinement_of(m2, hierarchies['modification']):
                self.modification_graph.add_edge(_mc_tuple(m2), _mc_tuple(m1))
            elif m2.refinement_of(m1, hierarchies['modification']):
                self.modification_graph.add_edge(_mc_tuple(m1), _mc_tuple(m2))

    def add_activity(self, activity_type):
        if activity_type not in self.activity_types:
            self.activity_types.append(activity_type)

    def add_active_state(self, activity_type, agent, evidence):
        agent_state = AgentState(agent, evidence)
        if activity_type in self.active_states:
            self.active_states[activity_type].append(agent_state)
        else:
            self.active_states[activity_type] = [agent_state]

    def add_inactive_state(self, activity_type, agent, evidence):
        agent_state = AgentState(agent, evidence)
        if activity_type in self.inactive_states:
            self.inactive_states[activity_type].append(agent_state)
        else:
            self.inactive_states[activity_type] = [agent_state]

    def get_active_forms(self):
        # TODO: handle activity types
        if self.active_states:
            states = []
            for k, v in self.active_states.items():
                states += v
            return states
        return None

    def get_inactive_forms(self):
        # TODO: handle activity types
        if self.inactive_states:
            states = []
            for k, v in self.inactive_states.items():
                states += v
            return states
        return None

    def add_modification(self, mc):
        mcc = ModCondition(mc.mod_type, mc.residue, mc.position, True)
        found = False
        for mod in self.modifications:
            if mcc.matches(mod):
                found = True
                break
        if not found:
            self.modifications.append(mcc)

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
    """A class representing Agent state without identifying a specific Agent.

    Attributes
    ----------
    bound_conditions : list[indra.statements.BoundCondition]
    mods : list[indra.statements.ModCondition]
    mutations : list[indra.statements.Mutation]
    location : indra.statements.location
    """
    def __init__(self, agent, evidence=None):
        self.bound_conditions = agent.bound_conditions
        self.mods = agent.mods
        self.mutations = agent.mutations
        self.location = agent.location
        self.evidence = evidence or []

    def apply_to(self, agent):
        """Apply this object's state to an Agent.

        Parameters
        ----------
        agent : indra.statements.Agent
            The agent to which the state should be applied
        """
        agent.bound_conditions = self.bound_conditions
        agent.mods = self.mods
        agent.mutations = self.mutations
        agent.location = self.location
        return self.evidence

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


def _get_statements_by_type(stmts, stmt_type):
    return [st for st in stmts if isinstance(st, stmt_type)]


def _get_graph_reductions(graph):
    """Return transitive reductions on a DAG.

    This is used to reduce the set of activities of a BaseAgent to the most
    specific one(s) possible. For instance, if a BaseAgent is know to have
    'activity', 'catalytic' and 'kinase' activity, then this function will
    return {'activity': 'kinase', 'catalytic': 'kinase', 'kinase': 'kinase'}
    as the set of reductions.
    """
    def frontier(g, nd):
        """Return the nodes after nd in the topological sort that are at the
        lowest possible level of the topological sort."""
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
    # This loop ensures that if a node n2 comes after node n1 in the topological
    # sort, and their frontiers are identical then n1 can be reduced to n2.
    # If their frontiers aren't identical, the reduction cannot be done.
    for i, n1 in enumerate(nodes_sort):
        for j, n2 in enumerate(nodes_sort):
            if i > j:
                continue
            if frontiers[i] == frontiers[j]:
                reductions[n1] = n2
    return reductions

def _mc_tuple(mc):
    return (mc.mod_type, mc.residue, mc.position)

