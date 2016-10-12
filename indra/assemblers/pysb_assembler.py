from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import itertools

from pysb import (Model, Monomer, Parameter, Rule, Annotation,
        ComponentDuplicateNameError, ComplexPattern, ReactionPattern, ANY,
        InvalidInitialConditionError)
from pysb.core import SelfExporter
import pysb.export

from indra import statements as ist
from indra.databases import context_client

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('pysb_assembler')

SelfExporter.do_export = False

# Here we define the types of INDRA statements that are meant to be
# assembled using the PySB assembler. If a type of statement appears
# in this list then we require that there is at least one default
# policy implemented to assemble that type of statement.
statement_whitelist = [ist.Phosphorylation, ist.Dephosphorylation,
                       ist.SelfModification, ist.Complex,
                       ist.Activation, ist.ActiveForm,
                       ist.RasGef, ist.RasGap, ist.Translocation]

def _n(name):
    """Return valid PySB name."""
    n = name.encode('ascii', errors='ignore').decode('ascii')
    n = re.sub('[^A-Za-z0-9_]', '_', n)
    n = re.sub(r'(^[0-9].*)', r'p\1', n)
    return n

def _is_whitelisted(stmt):
    """Return True if the statement type is in the whitelist."""
    for s in statement_whitelist:
        if isinstance(stmt, s):
            return True
    return False

# BaseAgent classes ####################################################

class _BaseAgentSet(object):
    """Container for a dict of BaseAgents with their names as keys."""
    def __init__(self):
        self.agents = {}

    def get_create_base_agent(self, agent):
        """Return base agent with given name, creating it if needed."""
        try:
            base_agent = self.agents[_n(agent.name)]
        except KeyError:
            base_agent = _BaseAgent(_n(agent.name))
            self.agents[_n(agent.name)] = base_agent

        # Handle bound conditions
        for bc in agent.bound_conditions:
            bound_base_agent = self.get_create_base_agent(bc.agent)
            bound_base_agent.create_site(get_binding_site_name(_n(agent.name)))
            base_agent.create_site(get_binding_site_name(_n(bc.agent.name)))

        # Handle modification conditions
        for mc in agent.mods:
            mod_site_name =\
                get_mod_site_name(mc.mod_type, mc.residue, mc.position)
            site_states = states[mc.mod_type]
            base_agent.create_site(mod_site_name, site_states)

        # Handle mutation conditions
        for mc in agent.mutations:
            if mc.residue_from is None:
                res_from = 'X'
            else:
                res_from = mc.residue_from
            mut_site_name = res_from + mc.position
            base_agent.create_site(mut_site_name, states=['WT'])
            if mc.residue_to is not None:
                base_agent.add_site_states(mut_site_name, [mc.residue_to])

        # Handle location condition
        if agent.location is not None:
            base_agent.create_site('loc', [agent.location])

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


class _BaseAgent(object):
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
        # The list of site/state configurations that lead to this agent
        # being active (where the agent is currently assumed to have only
        # one type of activity)
        self.active_forms = []
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
            self.active_forms.append(activity_pattern)
        else:
            self.inactive_forms.append(activity_pattern)

# Site/state information ###############################################

abbrevs = {
    'phosphorylation': 'phospho',
    'ubiquitination': 'ub',
    'farnesylation': 'farnesyl',
    'hydroxylation': 'hydroxyl',
    'acetylation': 'acetyl',
    'sumoylation': 'sumo',
    'glycosylation': 'glycosyl',
    'methylation': 'methyl',
    'modification': 'mod',
}

active_site_names = {
    'kinase': 'kin_site',
    'phosphatase': 'phos_site',
    'gtpbound': 'switch',
    'catalytic': 'cat_site',
    'transcription': 'trans_act',
    # For general molecular activity
    'activity': 'act'
}

states = {
    'phosphorylation': ['u', 'p'],
    'ubiquitination': ['n', 'y'],
    'farnesylation': ['n', 'y'],
    'hydroxylation': ['n', 'y'],
    'acetylation': ['n', 'y'],
    'sumoylation': ['n', 'y'],
    'glycosylation': ['n', 'y'],
    'methylation': ['n', 'y'],
    'modification': ['n', 'y'],
}

# The following dict specifies the default modification/binding site names for
# modifications resulting from a particular type of activity. For example, a
# protein with Kinase activity makes a modification of type "phospho" on its
# substrate, and a RasGTPase (with GtpBound activity) binds to a site of type
# "RBD" (Ras binding domain). This comes in handy for specifying
# Activation rules, where the modification site mediating the activation
# is not specified.
default_mod_site_names = {
    'kinase': 'phospho',
    'gtpbound': 'rbd',
    'phosphatase': 'phospho',
    'activity': 'act'
}


def get_binding_site_name(name):
    """Return a binding site name from a given agent name."""
    binding_site = _n(name).lower()
    return binding_site

def get_mod_site_name(mod_type, residue, position):
    """Return site names for a modification."""
    names = []
    if residue is None:
        mod_str = abbrevs[mod_type]
    else:
        mod_str = residue
    mod_pos = position if position is not None else ''
    name = ('%s%s' % (mod_str, mod_pos))
    return name

def get_active_forms(agent, agent_set):
    '''Returns all the patterns (dicts of site states) of an Agent
    that are known to be active.'''
    act_forms = agent_set[_n(agent.name)].active_forms
    if not act_forms:
        act_forms = [{}]
    return act_forms

def get_inactive_forms(agent, agent_set):
    '''Returns all the patterns (dicts of site states) of an Agent
    that are known to be inactive.'''
    inact_forms = agent_set[_n(agent.name)].inactive_forms
    if not inact_forms:
        inact_forms = [{}]
    return inact_forms

# PySB model elements ##################################################

def get_agent_rule_str(agent):
    """Construct a string from an Agent as part of a PySB rule name."""
    rule_str_list = [_n(agent.name)]
    for mod in agent.mods:
        mstr = abbrevs[mod.mod_type]
        if mod.residue is not None:
            mstr += mod.residue
        if mod.position is not None:
            mstr += mod.position
        rule_str_list.append('%s' % mstr)
    for mut in agent.mutations:
        if mut.residue_from is None:
            res_from = 'X'
        else:
            res_from = mut.residue_from
        mstr = res_from + mut.position
        if mut.residue_to is not None:
            mstr += mut.residue_to
        rule_str_list.append(mstr)
    if agent.bound_conditions:
        for b in agent.bound_conditions:
            if b.is_bound:
                rule_str_list.append(_n(b.agent.name))
            else:
                rule_str_list.append('n' + _n(b.agent.name))
    if agent.location is not None:
        rule_str_list.append(agent.location.replace(' ', '_'))
    rule_str = '_'.join(rule_str_list)
    return rule_str

def add_rule_to_model(model, rule):
    """Add a Rule to a PySB model and handle duplicate component errors."""
    try:
        model.add_component(rule)
    # If this rule is already in the model, issue a warning and continue
    except ComponentDuplicateNameError:
        msg = "Rule %s already in model! Skipping." % rule.name
        logger.warning(msg)


def get_create_parameter(model, name, value, unique=True):
    """Return parameter with given name, creating it if needed.

    If unique is false and the parameter exists, the value is not changed; if
    it does not exist, it will be created. If unique is true then upon conflict
    a number is added to the end of the parameter name.
    """
    parameter = model.parameters.get(name)

    if not unique and parameter is not None:
        return parameter

    if unique:
        pnum = 1
        while True:
            pname = name + '_%d' % pnum
            if model.parameters.get(pname) is None:
                break
            pnum += 1
    else:
        pname = name

    parameter = Parameter(pname, value)
    model.add_component(parameter)
    return parameter

def get_uncond_agent(agent):
    """Construct the unconditional state of an Agent.

    The unconditional Agent is a copy of the original agent but
    without any bound conditions and modification conditions.
    Mutation conditions, however, are preserved since they are static.
    """
    agent_uncond = ist.Agent(_n(agent.name), mutations=agent.mutations)
    return agent_uncond

def get_monomer_pattern(model, agent, extra_fields=None):
    """Construct a PySB MonomerPattern from an Agent."""
    pattern = get_site_pattern(agent)
    if extra_fields is not None:
        for k, v in extra_fields.items():
            pattern[k] = v

    # If a model is given, return the Monomer with the generated pattern,
    # otherwise just return the pattern
    monomer = model.monomers[_n(agent.name)]
    monomer_pattern = monomer(**pattern)
    return monomer_pattern

def get_site_pattern(agent):
    """Construct a dictionary of Monomer site states from an Agent.

    This crates the mapping to the associated PySB monomer from an
    INDRA Agent object."""
    pattern = {}
    # Handle bound conditions
    for bc in agent.bound_conditions:
        # Here we make the assumption that the binding site
        # is simply named after the binding partner
        if bc.is_bound:
            pattern[get_binding_site_name(_n(bc.agent.name))] = ANY
        else:
            pattern[get_binding_site_name(_n(bc.agent.name))] = None

    # Handle modifications
    for mod in agent.mods:
        mod_site_str = abbrevs[mod.mod_type]
        if mod.residue is not None:
            mod_site_str = mod.residue
        mod_pos_str = mod.position if mod.position is not None else ''
        mod_site = ('%s%s' % (mod_site_str, mod_pos_str))
        site_states = states[mod.mod_type]
        if mod.is_modified:
            pattern[mod_site] = site_states[1]
        else:
            pattern[mod_site] = site_states[0]

    # Handle mutations
    for mc in agent.mutations:
        if mc.residue_from is None:
            res_from = 'X'
        else:
            res_from = mc.residue_from
        mut_site_name = res_from + mc.position
        mut_site_state = mc.residue_to
        pattern[mut_site_name] = mut_site_state

    # Handle location
    if agent.location is not None:
        pattern['loc'] = agent.location

    return pattern


def set_base_initial_condition(model, monomer, value):
    """Set an initial condition for a monomer in its 'default' state."""
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

def set_extended_initial_condition(model, monomer, value=0):
    """Set an initial condition for monomers in "modified" state.

    This is useful when using downstream analysis that relies on reactions
    being active in the model. One example is BioNetGen-based reaction network
    diagram generation.
    """
    # Build up monomer pattern dict for default state
    sites_dict = {}
    for site in monomer.sites:
        if site in monomer.site_states:
            sites_dict[site] = monomer.site_states[site][-1]
        else:
            sites_dict[site] = None
    mp = monomer(**sites_dict)
    pname = monomer.name + '_0_mod'
    try:
        p = model.parameters[pname]
        p.value = value
    except KeyError:
        p = Parameter(pname, value)
        model.add_component(p)
        try:
            model.initial(mp, p)
        except InvalidInitialConditionError:
            pass

def get_annotation(component, db_name, db_ref):
    """Construct model Annotations for each component.

    Annotation formats follow guidelines at http://identifiers.org/.
    """
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
    elif db_name == 'IP':
        obj = url + 'interpro/%s' % db_ref
        pred = 'is'
    elif db_name == 'CHEBI':
        obj = url + 'chebi/CHEBI:%s' % db_ref
        pred = 'is'
    else:
        return None
    return Annotation(subj, obj, pred)

# PysbAssembler #######################################################

class UnknownPolicyError(Exception):
    pass

class PysbAssembler(object):
    """Assembler creating a PySB model from a set of INDRA Statements.

    Parameters
    ----------
    policies : Optional[Union[str, dict]]
        A string or dictionary that defines one or more assembly policies.

        If policies is a string, it defines a global assembly policy
        that applies to all Statement types.
        Example: contact_only, one_step

        A dictionary of policies has keys corresponding to Statement types
        and values to the policy to be applied to that type of Statement.
        For Statement types whose policy is undefined, the 'default'
        policy is applied.
        Example: {'phosphorylation': 'two_step'}

    Attributes
    ----------
    policies : dict
        A dictionary of policies that defines assembly policies for Statement
        types. It is assigned in the constructor.
    statements : list
        A list of INDRA statements to be assembled.
    model : pysb.Model
        A PySB model object that is assembled by this class.
    agent_set : _BaseAgentSet
        A set of BaseAgents used during the assembly process.
    """
    def __init__(self, policies=None):
        self.statements = []
        self.agent_set = None
        self.model = None
        if policies is None:
            self.policies = {'other': 'default'}
        elif isinstance(policies, basestring):
            self.policies = {'other': policies}
        else:
            self.policies = {'other': 'default'}
            self.policies.update(policies)

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        for stmt in stmts:
            if not self._statement_exists(stmt):
                self.statements.append(stmt)

    def make_model(self, policies=None, initial_conditions=True):
        """Assemble the PySB model from the collected INDRA Statements.

        This method assembles a PySB model from the set of INDRA Statements.
        The assembled model is both returned and set as the assembler's
        model argument.

        Parameters
        ----------
        policies : Optional[Union[str, dict]]
            A string or dictionary of policies, as defined in
            :py:class:`indra.assemblers.PysbAssembler`. This set of policies
            locally supersedes the default setting in the assembler. This
            is useful when this function is called multiple times with
            different policies.
        initial_conditions : Optional[bool]
            If True, default initial conditions are generated for the
            Monomers in the model.

        Returns
        -------
        model : pysb.Model
            The assembled PySB model object.
        """
        # Set local policies for this make_model call that overwrite
        # the global policies of the PySB assembler
        if policies is not None:
            global_policies = self.policies
            if isinstance(policies, str):
                local_policies = {'other': policies}
            else:
                local_policies = {'other': 'default'}
                local_policies.update(policies)
            self.policies = local_policies
        self.model = Model()
        self.agent_set = _BaseAgentSet()
        # Collect information about the monomers/self.agent_set from the
        # statements
        self._monomers()
        # Add the monomers to the model based on our BaseAgentSet
        for agent_name, agent in self.agent_set.items():
            m = Monomer(_n(agent_name), agent.sites, agent.site_states)
            self.model.add_component(m)
            for db_name, db_ref in agent.db_refs.items():
                a = get_annotation(m, db_name, db_ref)
                if a is not None:
                    self.model.add_annotation(a)
        # Iterate over the statements to generate rules
        self._assemble()
        # Add initial conditions
        if initial_conditions:
            self.add_default_initial_conditions()

        # If local policies were applied, revert to the global one
        if policies is not None:
            self.policies = global_policies

        return self.model

    def add_default_initial_conditions(self):
        """Set default initial conditions in the PySB model."""
        if self.model is None:
            return
        for m in self.model.monomers:
            set_base_initial_condition(self.model, m, 100.0)

    def set_context(self, cell_type):
        """Set protein expression data as initial conditions.

        This method uses :py:mod:`indra.databases.context_client` to get
        protein expression levels for a given cell type and set initial
        conditions for Monomers in the model accordingly.

        Parameters
        ----------
        cell_type : str
            Cell type name for which expression levels are queried.
            The cell type name follows the CCLE database conventions.

        Example: LOXIMVI_SKIN, BT20_BREAST
        """
        if self.model is None:
            return
        monomer_names = [m.name for m in self.model.monomers]
        res = context_client.get_protein_expression(monomer_names, cell_type)
        if not res:
            logger.warning('Could not get context for %s cell type.' %
                           cell_type)
            self.add_default_initial_conditions()
        monomers_found = []
        monomers_notfound = []
        for m in self.model.monomers:
            init = res.get(m.name)
            if init is not None:
                set_base_initial_condition(self.model, m, init[cell_type])
                monomers_found.append(m.name)
            else:
                set_base_initial_condition(self.model, m, 100.0)
                monomers_notfound.append(m.name)
        logger.info('Monomers set to %s context' % cell_type)
        logger.info('--------------------------------')
        for m in monomers_found:
            logger.info('%s' % m)
        if monomers_notfound:
            logger.info('')
            logger.info('Monomers not found in %s context' % cell_type)
            logger.info('-----------------------------------------')
            for m in monomers_notfound:
                logger.info('%s' % m)

    def print_model(self):
        """Print the assembled model as a PySB program string.

        This function is useful when the model needs to be passed as a string
        to another component.
        """
        model_str = pysb.export.export(self.model, 'pysb_flat')
        return model_str

    def save_model(self, file_name='pysb_model.py'):
        """Save the assembled model as a PySB program file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the model program code in.
            Default: pysb-model.py
        """
        if self.model is not None:
            model_str = self.print_model()
            with open(file_name, 'wt') as fh:
                fh.write(model_str)

    def save_rst(self, file_name='pysb_model.rst', module_name='pysb_module'):
        """Save the assembled model as an RST file for literate modeling.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the RST in.
            Default: pysb_model.rst
        module_name : Optional[str]
            The name of the python function defining the module.
            Default: pysb_module
        """
        if self.model is not None:
            with open(file_name, 'wt') as fh:
                fh.write('.. _%s:\n\n' % module_name)
                fh.write('Module\n======\n\n')
                fh.write('INDRA-assembled model\n---------------------\n\n')
                fh.write('::\n\n')
                model_str = pysb.export.export(self.model, 'pysb_flat')
                model_str = '\t' + model_str.replace('\n', '\n\t')
                fh.write(model_str)

    def _statement_exists(self, stmt):
        """Return True if the given Statement is in the assembler."""
        for s in self.statements:
            if stmt.matches(s):
                return True
        return False

    def _dispatch(self, stmt, stage, *args):
        """Construct and call an assembly function.

        This function constructs the name of the assembly function based on
        the type of statement, the corresponding policy and the stage
        of assembly. It then calls that function to perform the assembly
        task."""
        class_name = stmt.__class__.__name__
        try:
            policy = self.policies[class_name]
        except KeyError:
            policy = self.policies['other']
        func_name = '%s_%s_%s' % (class_name.lower(), stage, policy)
        func = globals().get(func_name)
        if func is None:
            # The specific policy is not implemented for the
            # given statement type.
            # We try to apply a default policy next.
            func_name = '%s_%s_default' % (class_name.lower(), stage)
            func = globals().get(func_name)
            if func is None:
                # The given statement type doesn't have a default
                # policy.
                raise UnknownPolicyError('%s function %s not defined' %
                                         (stage, func_name))
        return func(stmt, *args)

    def _monomers(self):
        """Calls the appropriate monomers method based on policies."""
        for stmt in self.statements:
            if _is_whitelisted(stmt):
                self._dispatch(stmt, 'monomers', self.agent_set)

    def _assemble(self):
        """Calls the appropriate assemble method based on policies."""
        for stmt in self.statements:
            if _is_whitelisted(stmt):
                self._dispatch(stmt, 'assemble', self.model, self.agent_set)


# COMPLEX ############################################################

def complex_monomers_one_step(stmt, agent_set):
    """In this (very simple) implementation, proteins in a complex are
    each given site names corresponding to each of the other members
    of the complex (lower case). So the resulting complex can be
    "fully connected" in that each member can be bound to
    all the others."""
    for i, member in enumerate(stmt.members):
        gene_mono = agent_set.get_create_base_agent(member)

        # Specify a binding site for each of the other complex members
        # bp = abbreviation for "binding partner"
        for j, bp in enumerate(stmt.members):
            # The protein doesn't bind to itstmt!
            if i == j:
                continue
            gene_mono.create_site(get_binding_site_name(bp.name))

complex_monomers_default = complex_monomers_one_step


def complex_assemble_one_step(stmt, model, agent_set):
    pairs = itertools.combinations(stmt.members, 2)
    for pair in pairs:
        agent1 = pair[0]
        agent2 = pair[1]
        param_name = agent1.name[0].lower() + \
                     agent2.name[0].lower() + '_bind'
        kf_bind = get_create_parameter(model, 'kf_' + param_name, 1e-6)
        kr_bind = get_create_parameter(model, 'kr_' + param_name, 1e-6)

        # Make a rule name
        rule_name = '_'.join([get_agent_rule_str(m) for m in pair])
        rule_name += '_bind'

        # Construct full patterns of each agent with conditions
        agent1_pattern = get_monomer_pattern(model, agent1)
        agent2_pattern = get_monomer_pattern(model, agent2)
        agent1_bs = get_binding_site_name(agent2.name)
        agent2_bs = get_binding_site_name(agent1.name)
        r = Rule(rule_name, agent1_pattern(**{agent1_bs: None}) + \
                            agent2_pattern(**{agent2_bs: None}) >>
                            agent1_pattern(**{agent1_bs: 1}) % \
                            agent2_pattern(**{agent2_bs: 1}),
                            kf_bind)
        add_rule_to_model(model, r)

        # In reverse reaction, assume that dissocition is unconditional

        agent1_uncond = get_uncond_agent(agent1)
        agent1_rule_str = get_agent_rule_str(agent1_uncond)
        monomer1_uncond = get_monomer_pattern(model, agent1_uncond)
        agent2_uncond = get_uncond_agent(agent2)
        agent2_rule_str = get_agent_rule_str(agent2_uncond)
        monomer2_uncond = get_monomer_pattern(model, agent2_uncond)
        rule_name = '%s_%s_dissociate' % (agent1_rule_str, agent2_rule_str)
        r = Rule(rule_name, monomer1_uncond(**{agent1_bs: 1}) % \
                            monomer2_uncond(**{agent2_bs: 1}) >>
                            monomer1_uncond(**{agent1_bs: None}) + \
                            monomer2_uncond(**{agent2_bs: None}),
                            kr_bind)
        add_rule_to_model(model, r)


def complex_assemble_multi_way(stmt, model, agent_set):
    # Get the rate parameter
    abbr_name = ''.join([m.name[0].lower() for m in stmt.members])
    kf_bind = get_create_parameter(model, 'kf_' + abbr_name + '_bind', 1e-6)
    kr_bind = get_create_parameter(model, 'kr_' + abbr_name + '_bind', 1e-6)

    # Make a rule name
    rule_name = '_'.join([get_agent_rule_str(m) for m in stmt.members])
    rule_name += '_bind'

    # Initialize the left and right-hand sides of the rule
    lhs = ReactionPattern([])
    rhs = ComplexPattern([], None)
    # We need a unique bond index for each pair of proteins in the
    # complex, resulting in n(n-1)/2 bond indices for a n-member complex.
    # We keep track of the bond indices using the bond_indices dict,
    # which maps each unique pair of members to a bond index.
    bond_indices = {}
    bond_counter = 1
    for i, member in enumerate(stmt.members):
        gene_name = member.name
        mono = model.monomers[gene_name]
        # Specify free and bound states for binding sites for each of
        # the other complex members
        # (bp = abbreviation for "binding partner")
        left_site_dict = {}
        right_site_dict = {}
        for j, bp in enumerate(stmt.members):
            bp_bs = get_binding_site_name(bp.name)
            # The protein doesn't bind to itstmt!
            if i == j:
                continue
            # Check to see if we've already created a bond index for these
            # two binding partners
            bp_set = frozenset([i, j])
            if bp_set in bond_indices:
                bond_ix = bond_indices[bp_set]
            # If we haven't see this pair of proteins yet, add a new bond
            # index to the dict
            else:
                bond_ix = bond_counter
                bond_indices[bp_set] = bond_ix
                bond_counter += 1
            # Fill in the entries for the site dicts
            left_site_dict[bp_bs] = None
            right_site_dict[bp_bs] = bond_ix

        # Add the pattern for the modifications of the member
        for mod in member.mods:
            if mod.residue is None:
                mod_str = abbrevs[mod.mod_type]
            else:
                mod_str = mod.residue
            mod_pos = mod.position if mod.position is not None else ''
            mod_site = ('%s%s' % (mod_str, mod_pos))
            left_site_dict[mod_site] = states[mod.mod_type][1]
            right_site_dict[mod_site] = states[mod.mod_type][1]

        # Add the pattern for the member being bound
        for bc in member.bound_conditions:
            bound_name = _n(bc.agent.name)
            bound_bs = get_binding_site_name(bound_name)
            gene_bs = get_binding_site_name(gene_name)
            if bc.is_bound:
                bound = model.monomers[bound_name]
                left_site_dict[bound_bs] = \
                    bond_counter
                right_site_dict[bound_bs] = \
                    bond_counter
                left_pattern = mono(**left_site_dict) % \
                                bound(**{gene_bs: bond_counter})
                right_pattern = mono(**right_site_dict) % \
                                bound(**{gene_bs: bond_counter})
                bond_counter += 1
            else:
                left_site_dict[bound_bs] = None
                right_site_dict[bound_bs] = None
                left_pattern = mono(**left_site_dict)
                right_pattern = mono(**right_site_dict)
        else:
            left_pattern = mono(**left_site_dict)
            right_pattern = mono(**right_site_dict)
        # Build up the left- and right-hand sides of the rule from
        # monomer patterns with the appropriate site dicts
        lhs = lhs + left_pattern
        rhs = rhs % right_pattern
    # Finally, create the rule and add it to the model
    rule_fwd = Rule(rule_name + '_fwd', lhs >> rhs, kf_bind)
    rule_rev = Rule(rule_name + '_rev', rhs >> lhs, kr_bind)
    add_rule_to_model(model, rule_fwd)
    add_rule_to_model(model, rule_rev)

complex_assemble_default = complex_assemble_one_step

# PHOSPHORYLATION ###################################################

def phosphorylation_monomers_interactions_only(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    enz.create_site(active_site_names['Kinase'])
    sub = agent_set.get_create_base_agent(stmt.sub)
    # See NOTE in monomers_one_step, below
    site_name = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(site_name, ('u', 'p'))


def phosphorylation_monomers_one_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    # NOTE: This assumes that a Phosphorylation statement will only ever
    # involve a single phosphorylation site on the substrate (typically
    # if there is more than one site, they will be parsed into separate
    # Phosphorylation statements, i.e., phosphorylation is assumed to be
    # distributive. If this is not the case, this assumption will need to
    # be revisited.
    site_name = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(site_name, ('u', 'p'))


def phosphorylation_monomers_two_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    site_name = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(site_name, ('u', 'p'))

    # Create site for binding the substrate
    enz.create_site(get_binding_site_name(sub.name))
    sub.create_site(get_binding_site_name(enz.name))

def phosphorylation_monomers_atp_dependent(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    site_name = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(site_name, ('u', 'p'))

    # Create site for binding the substrate
    enz.create_site(get_binding_site_name(sub.name))
    sub.create_site(get_binding_site_name(enz.name))

    # Make ATP base agent and create binding sites
    atp = agent_set.get_create_base_agent(ist.Agent('ATP'))
    atp.create_site('b')
    enz.create_site('ATP')


phosphorylation_monomers_default = phosphorylation_monomers_one_step


def phosphorylation_assemble_interactions_only(stmt, model, agent_set):
    if stmt.enz is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    kr_bind = get_create_parameter(model, 'kr_bind', 1.0, unique=False)

    enz = model.monomers[stmt.enz.name]
    sub = model.monomers[stmt.sub.name]

    # See NOTE in monomers_one_step
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)

    rule_name = '%s_phospho_%s_%s' % (rule_enz_str, rule_sub_str, site)
    active_site = active_site_names['Kinase']
    # Create a rule specifying that the substrate binds to the kinase at
    # its active site
    lhs = enz(**{active_site: None}) + sub(**{phos_site: None})
    rhs = enz(**{active_site: 1}) + sub(**{phos_site: 1})
    r_fwd = Rule(rule_name + '_fwd', lhs >> rhs, kf_bind)
    r_rev = Rule(rule_name + '_rev', rhs >> lhs, kr_bind)
    add_rule_to_model(model, r_fwd)
    add_rule_to_model(model, r_rev)


def phosphorylation_assemble_one_step(stmt, model, agent_set):
    if stmt.enz is None:
        return
    param_name = 'kf_' + stmt.enz.name[0].lower() + \
                    stmt.sub.name[0].lower() + '_phos'
    kf_phospho = get_create_parameter(model, param_name, 1e-6)

    # See NOTE in monomers_one_step
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    enz_pattern = get_monomer_pattern(model, stmt.enz)
    sub_unphos = get_monomer_pattern(model, stmt.sub,
        extra_fields={phos_site: 'u'})
    sub_phos = get_monomer_pattern(model, stmt.sub,
        extra_fields={phos_site: 'p'})

    enz_act_mods = get_active_forms(stmt.enz, agent_set)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    for i, am in enumerate(enz_act_mods):
        rule_name = '%s_phospho_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, phos_site, i + 1)
        r = Rule(rule_name,
                enz_pattern(am) + sub_unphos >>
                enz_pattern(am) + sub_phos,
                kf_phospho)
        add_rule_to_model(model, r)


def phosphorylation_assemble_two_step(stmt, model, agent_set):
    if stmt.enz is None:
        return
    sub_bs = get_binding_site_name(stmt.sub.name)
    enz_bound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: 1})
    enz_unbound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: None})
    sub_pattern = get_monomer_pattern(model, stmt.sub)

    param_name = ('kf_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kf_bind = get_create_parameter(model, param_name, 1e-6)
    param_name = ('kr_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kr_bind = get_create_parameter(model, param_name, 1e-3)
    param_name = ('kc_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_phos')
    kf_phospho = get_create_parameter(model, param_name, 1e-3)

    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    enz_act_mods = get_active_forms(stmt.enz, agent_set)
    enz_bs = get_binding_site_name(stmt.enz.name)
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    for i, am in enumerate(enz_act_mods):
        rule_name = '%s_phospho_bind_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, phos_site, i + 1)
        r = Rule(rule_name,
            enz_unbound(am) + \
            sub_pattern(**{phos_site: 'u', enz_bs: None}) >>
            enz_bound(am) % \
            sub_pattern(**{phos_site: 'u', enz_bs: 1}),
            kf_bind)
        add_rule_to_model(model, r)

        rule_name = '%s_phospho_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, phos_site, i + 1)
        r = Rule(rule_name,
            enz_bound(am) % \
                sub_pattern(**{phos_site: 'u', enz_bs: 1}) >>
            enz_unbound(am) + \
                sub_pattern(**{phos_site: 'p', enz_bs: None}),
            kf_phospho)
        add_rule_to_model(model, r)

    enz_uncond = get_uncond_agent(stmt.enz)
    enz_rule_str = get_agent_rule_str(enz_uncond)
    enz_mon_uncond = get_monomer_pattern(model, enz_uncond)
    sub_uncond = get_uncond_agent(stmt.sub)
    sub_rule_str = get_agent_rule_str(sub_uncond)
    sub_mon_uncond = get_monomer_pattern(model, sub_uncond)

    rule_name = '%s_dissoc_%s' % (enz_rule_str, sub_rule_str)
    r = Rule(rule_name, enz_mon_uncond(**{sub_bs: 1}) % \
             sub_mon_uncond(**{enz_bs: 1}) >>
             enz_mon_uncond(**{sub_bs: None}) + \
             sub_mon_uncond(**{enz_bs: None}), kr_bind)
    add_rule_to_model(model, r)

def phosphorylation_assemble_atp_dependent(stmt, model, agent_set):
    if stmt.enz is None:
        return
    # ATP
    atp = model.monomers['ATP']
    atp_bs = 'ATP'
    # ATP-bound enzyme
    enz_atp_bound = get_monomer_pattern(model, stmt.enz,
        extra_fields={atp_bs: 1})
    # ATP-free enzyme
    enz_atp_unbound = get_monomer_pattern(model, stmt.enz,
        extra_fields={atp_bs: None})
    # Substrate-bound enzyme
    sub_bs = get_binding_site_name(stmt.sub.name)
    enz_sub_bound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: 1})
    # Substrte and ATP-bound enzyme
    enz_sub_atp_bound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: 1, atp_bs: 2})
    enz_sub_atp_unbound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: None, atp_bs: None})
    # Substrate-free enzyme
    enz_sub_unbound = get_monomer_pattern(model, stmt.enz,
        extra_fields={sub_bs: None})
    # Enzyme active forms
    enz_act_forms = get_active_forms(stmt.enz, agent_set)
    enz_bs = get_binding_site_name(stmt.enz.name)
    # Unconditional enzyme
    enz_uncond = get_uncond_agent(stmt.enz)
    enz_rule_str = get_agent_rule_str(enz_uncond)
    enz_mon_uncond = get_monomer_pattern(model, enz_uncond)
    # Substrate
    sub_uncond = get_uncond_agent(stmt.sub)
    sub_rule_str = get_agent_rule_str(sub_uncond)
    sub_mon_uncond = get_monomer_pattern(model, sub_uncond)
    sub_pattern = get_monomer_pattern(model, stmt.sub)

    # Enzyme binding ATP
    param_name = ('kf_' + stmt.enz.name[0].lower() + '_atp_bind')
    kf_bind_atp = get_create_parameter(model, param_name, 1e-6)
    param_name = ('kr_' + stmt.enz.name[0].lower() + '_atp_bind')
    kr_bind_atp = get_create_parameter(model, param_name, 1e-6)
    for i, am in enumerate(enz_act_forms):
        rule_name = '%s_phospho_bind_atp_%d' % \
            (enz_rule_str, i + 1)
        r = Rule(rule_name,
            enz_atp_unbound(am) + atp(b=None) >>
            enz_atp_bound(am) %  atp(b=1), kf_bind_atp)
        add_rule_to_model(model, r)

    # Enzyme releasing ATP
    rule_name = '%s_phospho_dissoc_atp_%d' % \
        (enz_rule_str, i + 1)
    r = Rule(rule_name,
        enz_mon_uncond({atp_bs: 1}) % atp(b=1) >>
        enz_mon_uncond({atp_bs: None}) + atp(b=None), kr_bind_atp)
    add_rule_to_model(model, r)

    # Enzyme binding substrate
    param_name = ('kf_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kf_bind = get_create_parameter(model, param_name, 1e-6)
    param_name = ('kr_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kr_bind = get_create_parameter(model, param_name, 1e-3)
    param_name = ('kc_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_phos')
    kf_phospho = get_create_parameter(model, param_name, 1e-3)

    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    for i, am in enumerate(enz_act_forms):
        rule_name = '%s_phospho_bind_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, phos_site, i + 1)
        r = Rule(rule_name,
            enz_sub_unbound(am) + \
            sub_pattern(**{phos_site: 'u', enz_bs: None}) >>
            enz_sub_bound(am) % \
            sub_pattern(**{phos_site: 'u', enz_bs: 1}),
            kf_bind)
        add_rule_to_model(model, r)

    # Enzyme phosphorylating substrate
    for i, am in enumerate(enz_act_forms):
        rule_name = '%s_phospho_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, phos_site, i + 1)
        r = Rule(rule_name,
            enz_sub_atp_bound(am) % atp(b=2) % \
                sub_pattern(**{phos_site: 'u', enz_bs: 1}) >>
            enz_sub_atp_unbound(am) + atp(b=None) + \
                sub_pattern(**{phos_site: 'p', enz_bs: None}),
            kf_phospho)
        add_rule_to_model(model, r)

    # Enzyme dissodiating from substrate
    rule_name = '%s_dissoc_%s' % (enz_rule_str, sub_rule_str)
    r = Rule(rule_name, enz_mon_uncond(**{sub_bs: 1}) % \
             sub_mon_uncond(**{enz_bs: 1}) >>
             enz_mon_uncond(**{sub_bs: None}) + \
             sub_mon_uncond(**{enz_bs: None}), kr_bind)
    add_rule_to_model(model, r)

phosphorylation_assemble_default = phosphorylation_assemble_one_step

# CIS-AUTOPHOSPHORYLATION ###################################################

def autophosphorylation_monomers_interactions_only(stmt, agent_set):
    enz = agent_set.get_create_base_agent(stmt.enz)
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    enz.create_site(phos_site, ('u', 'p'))


def autophosphorylation_monomers_one_step(stmt, agent_set):
    enz = agent_set.get_create_base_agent(stmt.enz)
    # NOTE: This assumes that a Phosphorylation statement will only ever
    # involve a single phosphorylation site on the substrate (typically
    # if there is more than one site, they will be parsed into separate
    # Phosphorylation statements, i.e., phosphorylation is assumed to be
    # distributive. If this is not the case, this assumption will need to
    # be revisited.
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    enz.create_site(phos_site, ('u', 'p'))

autophosphorylation_monomers_default = autophosphorylation_monomers_one_step


def autophosphorylation_assemble_interactions_only(stmt, model, agent_set):
    stmt.assemble_one_step(model, agent_set)


def autophosphorylation_assemble_one_step(stmt, model, agent_set):
    param_name = 'kf_' + stmt.enz.name[0].lower() + '_autophos'
    kf_autophospho = get_create_parameter(model, param_name, 1e-3)

    # See NOTE in monomers_one_step
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    pattern_unphos = get_monomer_pattern(model, stmt.enz,
                                         extra_fields={phos_site: 'u'})
    pattern_phos = get_monomer_pattern(model, stmt.enz,
                                       extra_fields={phos_site: 'p'})
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_name = '%s_autophospho_%s_%s' % (rule_enz_str, rule_enz_str,
                                          phos_site)
    r = Rule(rule_name, pattern_unphos >> pattern_phos, kf_autophospho)
    add_rule_to_model(model, r)

autophosphorylation_assemble_default = autophosphorylation_assemble_one_step

# TRANSPHOSPHORYLATION ###################################################

def transphosphorylation_monomers_interactions_only(stmt, agent_set):
    enz = agent_set.get_create_base_agent(stmt.enz)
    # Assume there is exactly one bound_to species
    sub = agent_set.get_create_base_agent(stmt.enz)
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(phos_site, ('u', 'p'))


def transphosphorylation_monomers_one_step(stmt, agent_set):
    enz = agent_set.get_create_base_agent(stmt.enz)
    # NOTE: This assumes that a Phosphorylation statement will only ever
    # involve a single phosphorylation site on the substrate (typically
    # if there is more than one site, they will be parsed into separate
    # Phosphorylation statements, i.e., phosphorylation is assumed to be
    # distributive. If this is not the case, this assumption will need to
    # be revisited.
    sub = agent_set.get_create_base_agent(stmt.enz.bound_conditions[0].agent)
    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(phos_site, ('u', 'p'))

transphosphorylation_monomers_default = transphosphorylation_monomers_one_step


def transphosphorylation_assemble_interactions_only(stmt, model, agent_set):
    stmt.assemble_one_step(model, agent_set)


def transphosphorylation_assemble_one_step(stmt, model, agent_set):
    param_name = ('kf_' + stmt.enz.name[0].lower() +
                  _n(stmt.enz.bound_conditions[0].agent.name[0]).lower() +
                  '_transphos')
    kf = get_create_parameter(model, param_name, 1e-3)

    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    enz_pattern = get_monomer_pattern(model, stmt.enz)
    bound_agent = stmt.enz.bound_conditions[0].agent
    sub_unphos = get_monomer_pattern(model, bound_agent,
                                     extra_fields={phos_site: 'u'})
    sub_phos = get_monomer_pattern(model, bound_agent,
                                   extra_fields={phos_site: 'p'})

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_bound_str = get_agent_rule_str(bound_agent)
    rule_name = '%s_transphospho_%s_%s' % (rule_enz_str,
                                           rule_bound_str, phos_site)
    r = Rule(rule_name, enz_pattern % sub_unphos >> \
                    enz_pattern % sub_phos, kf)
    add_rule_to_model(model, r)

transphosphorylation_assemble_default = transphosphorylation_assemble_one_step

# ACTIVITYACTIVITY ######################################################

def activation_monomers_interactions_only(stmt, agent_set):
    subj = agent_set.get_create_base_agent(stmt.subj)
    obj = agent_set.get_create_base_agent(stmt.obj)
    if stmt.subj_activity is not None:
        subj.create_site(active_site_names[stmt.subj_activity])
        obj.create_site(active_site_names[stmt.obj_activity])
        obj.create_site(default_mod_site_names[stmt.subj_activity])
    else:
        subj.create_site(active_site_names[stmt.obj_activity])
        obj.create_site(active_site_names[stmt.obj_activity])
        obj.create_site(default_mod_site_names[stmt.obj_activity])

def activation_monomers_one_step(stmt, agent_set):
    subj = agent_set.get_create_base_agent(stmt.subj)
    obj = agent_set.get_create_base_agent(stmt.obj)
    if stmt.subj_activity is not None:
        subj.create_site(active_site_names[stmt.subj_activity],
                         ('inactive', 'active'))
    obj.create_site(active_site_names[stmt.obj_activity],
                    ('inactive', 'active'))

activation_monomers_default = activation_monomers_one_step


def activation_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    subj = model.monomers[stmt.subj.name]
    obj = model.monomers[stmt.obj.name]
    if stmt.subj_activity is not None:
        subj_active_site = active_site_names[stmt.subj_activity]
        obj_mod_site = default_mod_site_names[stmt.obj_activity]
    else:
        subj_active_site = active_site_names[stmt.obj_activity]
        obj_mod_site = default_mod_site_names[stmt.obj_activity]

    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_subj_str = get_agent_rule_str(stmt.subj)
    rule_name = '%s_%s_activates_%s_%s' %\
             (rule_subj_str, stmt.subj_activity, rule_obj_str,
              stmt.obj_activity)
    r = Rule(rule_name,
             subj(**{subj_active_site: None}) +
             obj(**{obj_mod_site: None}) >>
             subj(**{subj_active_site: 1}) %
             obj(**{obj_mod_site: 1}),
             kf_bind)
    add_rule_to_model(model, r)


def activation_assemble_one_step(stmt, model, agent_set):
    if stmt.subj_activity is not None:
        subj_pattern = get_monomer_pattern(model, stmt.subj,
            extra_fields={active_site_names[stmt.subj_activity]: 'active'})
    else:
        subj_pattern = get_monomer_pattern(model, stmt.subj)

    obj_inactive = get_monomer_pattern(model, stmt.obj,
        extra_fields={active_site_names[stmt.obj_activity]: 'inactive'})
    obj_active = get_monomer_pattern(model, stmt.obj,
        extra_fields={active_site_names[stmt.obj_activity]: 'active'})

    param_name = 'kf_' + stmt.subj.name[0].lower() + \
                        stmt.obj.name[0].lower() + '_act'
    kf_one_step_activate = \
                   get_create_parameter(model, param_name, 1e-6)

    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_subj_str = get_agent_rule_str(stmt.subj)
    rule_name = '%s_%s_activates_%s_%s' % \
        (rule_subj_str, stmt.subj_activity, rule_obj_str,
         stmt.obj_activity)

    if stmt.is_activation:
        r = Rule(rule_name,
            subj_pattern + obj_inactive >> subj_pattern + obj_active,
            kf_one_step_activate)
    else:
        r = Rule(rule_name,
            subj_pattern + obj_active >> subj_pattern + obj_inactive,
            kf_one_step_activate)

    add_rule_to_model(model, r)

activation_assemble_default = activation_assemble_one_step

# DEPHOSPHORYLATION #####################################################

def dephosphorylation_monomers_interactions_only(stmt, agent_set):
    if stmt.enz is None:
        return
    phos = agent_set.get_create_base_agent(stmt.enz)
    phos.create_site(active_site_names['phosphatase'])
    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub = agent_set.get_create_base_agent(stmt.sub)
    sub.create_site(dephos_site, ('u', 'p'))


def dephosphorylation_monomers_one_step(stmt, agent_set):
    if stmt.enz is None:
        return
    phos = agent_set.get_create_base_agent(stmt.enz)
    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub = agent_set.get_create_base_agent(stmt.sub)
    sub.create_site(dephos_site, ('u', 'p'))


def dephosphorylation_monomers_two_step(stmt, agent_set):
    if stmt.enz is None:
        return
    phos = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    sub.create_site(dephos_site, ('u', 'p'))

    # Create site for binding the substrate
    phos.create_site(get_binding_site_name(sub.name))
    sub.create_site(get_binding_site_name(phos.name))

dephosphorylation_monomers_default = dephosphorylation_monomers_one_step


def dephosphorylation_assemble_interactions_only(stmt, model, agent_set):
    if stmt.enz is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    phos = model.monomers[stmt.enz.name]
    sub = model.monomers[stmt.sub.name]
    phos_site = active_site_names['phosphatase']
    # See NOTE in Phosphorylation.monomers_one_step
    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    r = Rule('%s_dephospho_%s_%s' %
             (rule_enz_str, rule_sub_str, phos_site),
             phos(**{phos_site: None}) + sub(**{dephos_site: None}) >>
             phos(**{phos_site: 1}) + sub(**{dephos_site: 1}),
             kf_bind)
    add_rule_to_model(model, r)


def dephosphorylation_assemble_one_step(stmt, model, agent_set):
    if stmt.enz is None:
        return
    param_name = 'kf_' + stmt.enz.name[0].lower() + \
                stmt.sub.name[0].lower() + '_dephos'
    kf_dephospho = get_create_parameter(model, param_name, 1e-6)

    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)
    phos_pattern = get_monomer_pattern(model, stmt.enz)
    sub_phos = get_monomer_pattern(model, stmt.sub,
                                   extra_fields={dephos_site: 'p'})
    sub_unphos = get_monomer_pattern(model, stmt.sub,
                                     extra_fields={dephos_site: 'u'})

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    r = Rule('%s_dephospho_%s_%s' %
             (rule_enz_str, rule_sub_str, dephos_site),
             phos_pattern + sub_phos >>
             phos_pattern + sub_unphos,
             kf_dephospho)
    add_rule_to_model(model, r)


def dephosphorylation_assemble_two_step(stmt, model, agent_set):
    if stmt.enz is None:
        return
    sub_bs = get_binding_site_name(stmt.sub.name)
    phos_bs = get_binding_site_name(stmt.enz.name)
    phos_bound = get_monomer_pattern(model, stmt.enz,
                                     extra_fields={sub_bs: 1})
    phos_unbound = get_monomer_pattern(model, stmt.enz,
                                       extra_fields={sub_bs: None})
    sub_pattern = get_monomer_pattern(model, stmt.sub)

    param_name = 'kf_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_bind'
    kf_bind = get_create_parameter(model, param_name, 1e-6)
    param_name = 'kr_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_bind'
    kr_bind = get_create_parameter(model, param_name, 1e-3)
    param_name = 'kc_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_dephos'
    kf_phospho = get_create_parameter(model, param_name, 1e-3)

    dephos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    phos_act_mods = get_active_forms(stmt.enz, agent_set)
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    for i, am in enumerate(phos_act_mods):
        rule_name = '%s_dephos_bind_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, dephos_site, i + 1)
        r = Rule(rule_name,
            phos_unbound(am) + \
            sub_pattern(**{dephos_site: 'p', phos_bs: None}) >>
            phos_bound(am) % \
            sub_pattern(**{dephos_site: 'p', phos_bs: 1}),
            kf_bind, kr_bind)
        add_rule_to_model(model, r)

        rule_name = '%s_dephos_%s_%s_%d' % \
            (rule_enz_str, rule_sub_str, dephos_site, i + 1)
        r = Rule(rule_name,
            phos_bound(am) % \
                sub_pattern(**{dephos_site: 'p', phos_bs: 1}) >>
            phos_unbound(am) + \
                sub_pattern(**{dephos_site: 'u', phos_bs: None}),
            kf_phospho)
        add_rule_to_model(model, r)

    enz_uncond = get_uncond_agent(stmt.enz)
    enz_rule_str = get_agent_rule_str(enz_uncond)
    enz_mon_uncond = get_monomer_pattern(model, enz_uncond)
    sub_uncond = get_uncond_agent(stmt.sub)
    sub_rule_str = get_agent_rule_str(sub_uncond)
    sub_mon_uncond = get_monomer_pattern(model, sub_uncond)

    rule_name = '%s_dissoc_%s' % (enz_rule_str, sub_rule_str)
    r = Rule(rule_name, enz_mon_uncond(**{sub_bs: 1}) % \
             sub_mon_uncond(**{phos_bs: 1}) >>
             enz_mon_uncond(**{sub_bs: None}) + \
             sub_mon_uncond(**{phos_bs: None}), kr_bind)
    add_rule_to_model(model, r)

dephosphorylation_assemble_default = dephosphorylation_assemble_one_step

# RASGEF #####################################################

def rasgef_monomers_interactions_only(stmt, agent_set):
    gef = agent_set.get_create_base_agent(stmt.gef)
    gef.create_site('gef_site')
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('p_loop')


def rasgef_monomers_one_step(stmt, agent_set):
    gef = agent_set.get_create_base_agent(stmt.gef)
    gef.create_site(stmt.gef_activity, ('inactive', 'active'))
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('GtpBound', ('inactive', 'active'))

rasgef_monomers_default = rasgef_monomers_one_step


def rasgef_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    gef = model.monomers[stmt.gef.name]
    ras = model.monomers[stmt.ras.name]
    rule_gef_str = get_agent_rule_str(stmt.gef)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_activates_%s' %
             (rule_gef_str, rule_ras_str),
             gef(**{'gef_site': None}) +
             ras(**{'p_loop': None}) >>
             gef(**{'gef_site': 1}) +
             ras(**{'p_loop': 1}),
             kf_bind)
    add_rule_to_model(model, r)


def rasgef_assemble_one_step(stmt, model, agent_set):
    gef_pattern = get_monomer_pattern(model, stmt.gef,
        extra_fields={stmt.gef_activity: 'active'})
    ras_inactive = get_monomer_pattern(model, stmt.ras,
        extra_fields={'GtpBound': 'inactive'})
    ras_active = get_monomer_pattern(model, stmt.ras,
        extra_fields={'GtpBound': 'active'})

    param_name = 'kf_' + stmt.gef.name[0].lower() + \
                    stmt.ras.name[0].lower() + '_gef'
    kf_gef = get_create_parameter(model, param_name, 1e-6)

    rule_gef_str = get_agent_rule_str(stmt.gef)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_activates_%s' %
             (rule_gef_str, rule_ras_str),
             gef_pattern + ras_inactive >>
             gef_pattern + ras_active,
             kf_gef)
    add_rule_to_model(model, r)

rasgef_assemble_default = rasgef_assemble_one_step

# RASGAP ####################################################

def rasgap_monomers_interactions_only(stmt, agent_set):
    gap = agent_set.get_create_base_agent(stmt.gap)
    gap.create_site('gap_site')
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('gtp_site')


def rasgap_monomers_one_step(stmt, agent_set):
    gap = agent_set.get_create_base_agent(stmt.gap)
    gap.create_site(stmt.gap_activity, ('inactive', 'active'))
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('GtpBound', ('inactive', 'active'))

rasgap_monomers_default = rasgap_monomers_one_step


def rasgap_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    gap = model.monomers[stmt.gap.name]
    ras = model.monomers[stmt.ras.name]
    rule_gap_str = get_agent_rule_str(stmt.gap)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_inactivates_%s' %
             (rule_gap_str, rule_ras_str),
             gap(**{'gap_site': None}) +
             ras(**{'gtp_site': None}) >>
             gap(**{'gap_site': 1}) +
             ras(**{'gtp_site': 1}),
             kf_bind)
    add_rule_to_model(model, r)


def rasgap_assemble_one_step(stmt, model, agent_set):
    gap_pattern = get_monomer_pattern(model, stmt.gap,
        extra_fields={stmt.gap_activity: 'active'})
    ras_inactive = get_monomer_pattern(model, stmt.ras,
        extra_fields={'GtpBound': 'inactive'})
    ras_active = get_monomer_pattern(model, stmt.ras,
        extra_fields={'GtpBound': 'active'})

    param_name = 'kf_' + stmt.gap.name[0].lower() + \
                    stmt.ras.name[0].lower() + '_gap'
    kf_gap = get_create_parameter(model, param_name, 1e-6)

    rule_gap_str = get_agent_rule_str(stmt.gap)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_deactivates_%s' %
             (rule_gap_str, rule_ras_str),
             gap_pattern + ras_active >>
             gap_pattern + ras_inactive,
             kf_gap)
    add_rule_to_model(model, r)

rasgap_assemble_default = rasgap_assemble_one_step

# ACTIVEFORM ############################################

def activeform_monomers_interactions_only(stmt, agent_set):
    pass


def activeform_monomers_one_step(stmt, agent_set):
    agent = agent_set.get_create_base_agent(stmt.agent)
    site_conditions = get_site_pattern(stmt.agent)

    # Add this activity pattern explicitly to the agent's list
    # of active states
    agent.add_activity_form(site_conditions, stmt.is_active)

activeform_monomers_default = activeform_monomers_one_step


def activeform_assemble_interactions_only(stmt, model, agent_set):
    pass


def activeform_assemble_one_step(stmt, model, agent_set):
    pass

activeform_assemble_default = activeform_assemble_one_step

# RASGTPACTIVITIACTIVITY ######################################
def rasgtpactivation_monomers_default(stmt, agent_set):
    pass

def rasgtpactivation_assemble_default(stmt, model, agent_set):
    pass

# TRANSLOCATION ###############################################
def translocation_monomers_default(stmt, agent_set):
    # Skip if either from or to locations are missing
    if stmt.from_location is None or stmt.to_location is None:
        return
    agent = agent_set.get_create_base_agent(stmt.agent)
    agent.create_site('loc', [stmt.from_location, stmt.to_location])

def translocation_assemble_default(stmt, model, agent_set):
    if stmt.from_location is None or stmt.to_location is None:
        return
    param_name = 'kf_%s_%s_%s' % (_n(stmt.agent.name).lower(),
                                  stmt.from_location, stmt.to_location)
    kf_trans = get_create_parameter(model, param_name, 1.0, unique=True)
    monomer = model.monomers[_n(stmt.agent.name)]
    rule_agent_str = get_agent_rule_str(stmt.agent)
    rule_name = '%s_translocates_%s_to_%s' % (rule_agent_str,
                                              stmt.from_location,
                                              stmt.to_location)
    agent_from = get_monomer_pattern(model, stmt.agent,
                                     extra_fields={'loc': stmt.from_location})
    agent_to = get_monomer_pattern(model, stmt.agent,
                                   extra_fields={'loc': stmt.to_location})
    r = Rule(rule_name, agent_from >> agent_to, kf_trans)
    add_rule_to_model(model, r)


