import logging
import itertools
from indra.statements import *
from indra.util import fast_deepcopy
from .base_agents import BaseAgentSet

logger = logging.getLogger(__name__)


class PysbPreassembler(object):
    """Pre-assemble Statements in preparation for PySB assembly.

    Parameters
    ----------
    stmts : list[indra.statements.Statement]
        A list of Statements to assemble
    """
    def __init__(self, stmts=None):
        if not stmts:
            stmts = []
        self.statements = stmts
        self.agent_set = BaseAgentSet()

    def add_statements(self, stmts):
        """Add a list of Statements for assembly."""
        self.statements += stmts

    def _gather_active_forms(self):
        """Collect all the active forms of each Agent in the Statements."""
        for stmt in self.statements:
            if isinstance(stmt, ActiveForm):
                base_agent = self.agent_set.get_create_base_agent(stmt.agent)
                # Handle the case where an activity flag is set
                agent_to_add = stmt.agent
                if stmt.agent.activity:
                    new_agent = fast_deepcopy(stmt.agent)
                    new_agent.activity = None
                    agent_to_add = new_agent
                base_agent.add_activity_form(agent_to_add, stmt.is_active)

    def replace_activities(self):
        """Replace ative flags with Agent states when possible."""
        logger.debug('Running PySB Preassembler replace activities')
        # TODO: handle activity hierarchies
        new_stmts = []

        def has_agent_activity(stmt):
            """Return True if any agents in the Statement have activity."""
            for agent in stmt.agent_list():
                if isinstance(agent, Agent) and agent.activity is not None:
                    return True
            return False
        # First collect all explicit active forms
        self._gather_active_forms()
        # Iterate over all statements
        for j, stmt in enumerate(self.statements):
            logger.debug('%d/%d %s' % (j + 1, len(self.statements), stmt))
            # If the Statement doesn't have any activities, we can just
            # keep it and move on
            if not has_agent_activity(stmt):
                new_stmts.append(stmt)
                continue
            stmt_agents = stmt.agent_list()
            num_agents = len(stmt_agents)
            # Make a list with an empty list for each Agent so that later
            # we can build combinations of Agent forms
            agent_forms = [[] for a in stmt_agents]
            for i, agent in enumerate(stmt_agents):
                # This is the case where there is an activity flag on an
                # Agent which we will attempt to replace with an explicit
                # active form
                if agent is not None and isinstance(agent, Agent) and \
                        agent.activity is not None:
                    base_agent = self.agent_set.get_create_base_agent(agent)
                    # If it is an "active" state
                    if agent.activity.is_active:
                        active_forms = base_agent.active_forms
                        # If no explicit active forms are known then we use
                        # the generic one
                        if not active_forms:
                            active_forms = [agent]
                    # If it is an "inactive" state
                    else:
                        active_forms = base_agent.inactive_forms
                        # If no explicit inactive forms are known then we use
                        # the generic one
                        if not active_forms:
                            active_forms = [agent]
                    # We now iterate over the active agent forms and create
                    # new agents
                    for af in active_forms:
                        new_agent = fast_deepcopy(agent)
                        self._set_agent_context(af, new_agent)
                        agent_forms[i].append(new_agent)
                # Otherwise we just copy over the agent as is
                else:
                    agent_forms[i].append(agent)
            # Now create all possible combinations of the agents and create new
            # statements as needed
            agent_combs = itertools.product(*agent_forms)
            for agent_comb in agent_combs:
                new_stmt = fast_deepcopy(stmt)
                new_stmt.set_agent_list(agent_comb)
                new_stmts.append(new_stmt)
        self.statements = new_stmts

    def add_reverse_effects(self):
        """Add Statements for the reverse effects of some Statements.

        For instance, if a protein is phosphorylated but never dephosphorylated
        in the model, we add a generic dephosphorylation here. This step is
        usually optional in the assembly process.
        """
        # TODO: generalize to other modification sites
        pos_mod_sites = {}
        neg_mod_sites = {}
        syntheses = []
        degradations = []
        for stmt in self.statements:
            if isinstance(stmt, Phosphorylation):
                agent = stmt.sub.name
                try:
                    pos_mod_sites[agent].append((stmt.residue, stmt.position))
                except KeyError:
                    pos_mod_sites[agent] = [(stmt.residue, stmt.position)]
            elif isinstance(stmt, Dephosphorylation):
                agent = stmt.sub.name
                try:
                    neg_mod_sites[agent].append((stmt.residue, stmt.position))
                except KeyError:
                    neg_mod_sites[agent] = [(stmt.residue, stmt.position)]
            elif isinstance(stmt, Influence):
                if stmt.overall_polarity() == 1:
                    syntheses.append(stmt.obj.name)
                elif stmt.overall_polarity() == -1:
                    degradations.append(stmt.obj.name)
            elif isinstance(stmt, IncreaseAmount):
                syntheses.append(stmt.obj.name)
            elif isinstance(stmt, DecreaseAmount):
                degradations.append(stmt.obj.name)

        new_stmts = []
        for agent_name, pos_sites in pos_mod_sites.items():
            neg_sites = neg_mod_sites.get(agent_name, [])
            no_neg_site = set(pos_sites).difference(set(neg_sites))
            for residue, position in no_neg_site:
                st = Dephosphorylation(Agent('phosphatase'),
                                           Agent(agent_name),
                                           residue, position)
                new_stmts.append(st)
        for agent_name in syntheses:
            if agent_name not in degradations:
                st = DecreaseAmount(None, Agent(agent_name))
                new_stmts.append(st)

        self.statements += new_stmts

    @staticmethod
    def _set_agent_context(from_agent, to_agent):
        if not isinstance(from_agent, Agent) or \
                not isinstance(to_agent, Agent):
            return
        def add_no_duplicate(from_lst, to_lst):
            for fm in from_lst:
                found = False
                for tm in to_lst:
                    if fm.matches(tm):
                        found = True
                        break
                if not found:
                    to_lst.append(fm)
            return to_lst
        # TODO: what can we do about semantic conflicts here like the same
        # bound condition with True/False is_bound appearing in the
        # two contexts?
        to_agent.bound_conditions = \
            add_no_duplicate(to_agent.bound_conditions,
                             from_agent.bound_conditions)
        to_agent.mods = add_no_duplicate(to_agent.mods, from_agent.mods)
        to_agent.mutations = add_no_duplicate(to_agent.mutations,
                                              from_agent.mutations)
        to_agent.location = from_agent.location
        to_agent.activity = from_agent.activity

