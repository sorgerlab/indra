from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import itertools
from indra.util import fast_deepcopy

from pysb import (Model, Monomer, Parameter, Expression, Observable, Rule,
        Annotation, ComponentDuplicateNameError, ComplexPattern,
        ReactionPattern, ANY, WILD, InvalidInitialConditionError)
from pysb.core import SelfExporter
import pysb.export

from indra import statements as ist
from indra.databases import context_client, get_identifiers_url
from indra.preassembler.hierarchy_manager import entity_hierarchy as enth
from indra.tools.expand_families import _agent_from_uri

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
statement_whitelist = [ist.Modification, ist.SelfModification, ist.Complex,
                       ist.RegulateActivity, ist.ActiveForm,
                       ist.Gef, ist.Gap, ist.Translocation,
                       ist.IncreaseAmount, ist.DecreaseAmount,
                       ist.Conversion]

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
        site_name = get_mod_site_name(mc.mod_type,
                                      mc.residue, mc.position)
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
    'ribosylation': 'ribosyl',
    'geranylgeranylation': 'geranylgeranyl',
    'palmitoylation': 'palmitoyl',
    'myristoylation': 'myristoyl',
    'modification': 'mod',
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
    'geranylgeranylation': ['n', 'y'],
    'palmitoylation': ['n', 'y'],
    'myristoylation': ['n', 'y'],
    'ribosylation': ['n', 'y'],
    'modification': ['n', 'y'],
}

mod_acttype_map = {
    ist.Phosphorylation: 'kinase',
    ist.Dephosphorylation: 'phosphatase',
    ist.Hydroxylation: 'catalytic',
    ist.Dehydroxylation: 'catalytic',
    ist.Sumoylation: 'catalytic',
    ist.Desumoylation: 'catalytic',
    ist.Acetylation: 'catalytic',
    ist.Deacetylation: 'catalytic',
    ist.Glycosylation: 'catalytic',
    ist.Deglycosylation: 'catalytic',
    ist.Ribosylation: 'catalytic',
    ist.Deribosylation: 'catalytic',
    ist.Ubiquitination: 'catalytic',
    ist.Deubiquitination: 'catalytic',
    ist.Farnesylation: 'catalytic',
    ist.Defarnesylation: 'catalytic',
    ist.Palmitoylation: 'catalytic',
    ist.Depalmitoylation: 'catalytic',
    ist.Myristoylation: 'catalytic',
    ist.Demyristoylation: 'catalytic',
    ist.Geranylgeranylation: 'catalytic',
    ist.Degeranylgeranylation: 'catalytic',
    ist.Methylation: 'catalytic',
    ist.Demethylation: 'catalytic',
}


def get_binding_site_name(agent):
    """Return a binding site name from a given agent."""
    # Try to construct a binding site name based on parent
    grounding = agent.get_grounding()
    if grounding != (None, None):
        uri = enth.get_uri(grounding[0], grounding[1])
        # Get highest level parents in hierarchy
        parents = enth.get_parents(uri, 'top')
        if parents:
            # Choose the first parent if there are more than one
            parent_uri = sorted(list(parents))[0]
            parent_agent = _agent_from_uri(parent_uri)
            binding_site = _n(parent_agent.name).lower()
            return binding_site
    # Fall back to Agent's own name if one from parent can't be constructed
    binding_site = _n(agent.name).lower()
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
        res_from = mut.residue_from if mut.residue_from else 'mut'
        res_to = mut.residue_to if mut.residue_to else 'X'
        if mut.position is None:
            mut_site_name = res_from
        else:
            mut_site_name = res_from + mut.position
        mstr = mut_site_name + res_to
        rule_str_list.append(mstr)
    if agent.bound_conditions:
        for b in agent.bound_conditions:
            if b.is_bound:
                rule_str_list.append(_n(b.agent.name))
            else:
                rule_str_list.append('n' + _n(b.agent.name))
    if agent.location is not None:
        rule_str_list.append(_n(agent.location))
    rule_str = '_'.join(rule_str_list)
    return rule_str


def add_rule_to_model(model, rule, annotations=None):
    """Add a Rule to a PySB model and handle duplicate component errors."""
    try:
        model.add_component(rule)
        # If the rule was actually added, also add the annotations
        if annotations:
            model.annotations += annotations
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
    norm_name = _n(name)
    parameter = model.parameters.get(norm_name)

    if not unique and parameter is not None:
        return parameter

    if unique:
        pnum = 1
        while True:
            pname = norm_name + '_%d' % pnum
            if model.parameters.get(pname) is None:
                break
            pnum += 1
    else:
        pname = norm_name

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


def grounded_monomer_patterns(model, agent):
    """Get monomer patterns for the agent accounting for grounding information.
    """
    # Iterate over all model annotations to identify the monomer associated
    # with this agent
    monomer = None
    for ann in model.annotations:
        if monomer:
            break
        if not ann.predicate == 'is':
            continue
        if not isinstance(ann.subject, Monomer):
            continue
        (ns, id) = parse_identifiers_url(ann.object)
        if ns is None and id is None:
            continue
        # We now have an identifiers.org namespace/ID for a given monomer;
        # we check to see if there is a matching identifier in the db_refs
        # for this agent
        for db_ns, db_id in agent.db_refs.items():
            # We've found a match! Return first match
            # FIXME Could also update this to check for alternative
            # FIXME matches, or make sure that all grounding IDs match,
            # FIXME etc.
            if db_ns == ns and db_id == id:
                monomer = ann.subject
                break
    # We looked at all the annotations in the model and didn't find a
    # match
    if monomer is None:
        logger.info('No monomer found corresponding to agent %s' % agent)
        return
    # Now that we have a monomer for the agent, look for site/state
    # combinations corresponding to the state of the agent.  For every one of
    # the modifications specified in the agent signature, check to see if it
    # can be satisfied based on the agent's annotations.  For every one we find
    # that is consistent, we yield it--there may be more than one.
    # FIXME
    # Create a list of tuples, each one representing the site conditions
    # that can satisfy a particular agent condition. Each entry in the list
    # will contain a list of dicts associated with a particular mod/activity
    # condition. Each dict will represent a site/state combination satisfying
    # the constraints imposed by that mod/activity condition.
    sc_list = []
    for mod in agent.mods:
        # Find all site/state combinations that have the appropriate
        # modification type
        # As we iterate, build up a dict identifying the annotations of
        # particular sites
        mod_sites = {}
        res_sites = set([])
        pos_sites = set([])
        for ann in monomer.site_annotations:
            # Don't forget to handle Nones!
            if ann.predicate == 'is_modification' and \
               ann.object == mod.mod_type:
                site_state = ann.subject
                assert isinstance(site_state, tuple)
                assert len(site_state) == 2
                mod_sites[site_state[0]] = site_state[1]
            elif ann.predicate == 'is_residue' and \
                 ann.object == mod.residue:
                res_sites.add(ann.subject)
            elif ann.predicate == 'is_position' and \
                 ann.object == mod.position:
                pos_sites.add(ann.subject)
        # If the residue field of the agent is specified,
        viable_sites = set(mod_sites.keys())
        if mod.residue is not None:
            viable_sites = viable_sites.intersection(res_sites)
        if mod.position is not None:
            viable_sites = viable_sites.intersection(pos_sites)
        # If there are no viable sites annotated in the model matching the
        # available info in the mod condition, then we won't be able to
        # satisfy the conditions on this agent
        if not viable_sites:
            return
        # Otherwise, update the 
        # If there are any sites left after we subject them to residue
        # and position constraints, then return the relevant monomer patterns!
        pattern_list = []
        for site_name in viable_sites:
            pattern_list.append({site_name: (mod_sites[site_name], WILD)})
        sc_list.append(pattern_list)
    # Now check for monomer patterns satisfying the agent's activity condition
    if agent.activity:
        # Iterate through annotations with this monomer as the subject
        # and a has_active_pattern or has_inactive_pattern relationship
        # FIXME: Currently activity type is not annotated/checked
        # FIXME act_type = agent.activity.activity_type
        rel_type = 'has_active_pattern' if agent.activity.is_active \
                                        else 'has_inactive_pattern'
        active_form_list = []
        for ann in model.annotations:
            if ann.subject == monomer and ann.predicate == rel_type:
                # The annotation object contains the active/inactive pattern
                active_form_list.append(ann.object)
        sc_list.append(active_form_list)
    # Now that we've got a list of conditions
    for pattern_combo in itertools.product(*sc_list):
        mp_sc = {}
        for pattern in pattern_combo:
            mp_sc.update(pattern)
        if mp_sc:
            yield monomer(**mp_sc)
    if not sc_list:
        yield monomer()

def rules_with_annotation(model, monomer_name, predicate):
    rules = []
    for ann in model.annotations:
        if not ann.predicate == predicate:
            continue
        if ann.object == monomer_name:
            rules.append(model.rules[ann.subject])
    return rules


def get_monomer_pattern(model, agent, extra_fields=None):
    """Construct a PySB MonomerPattern from an Agent."""
    try:
        monomer = model.monomers[_n(agent.name)]
    except KeyError as e:
        logger.warning('Monomer with name %s not found in model' %
                       _n(agent.name))
        return None
    # Get the agent site pattern
    pattern = get_site_pattern(agent)
    if extra_fields is not None:
        for k, v in extra_fields.items():
            pattern[k] = v
    # If a model is given, return the Monomer with the generated pattern,
    # otherwise just return the pattern
    try:
        monomer_pattern = monomer(**pattern)
    except Exception as e:
        logger.info("Invalid site pattern %s for monomer %s" %
                      (pattern, monomer))
        return None
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
            pattern[get_binding_site_name(bc.agent)] = ANY
        else:
            pattern[get_binding_site_name(bc.agent)] = None

    # Handle modifications
    for mod in agent.mods:
        mod_site_str = abbrevs[mod.mod_type]
        if mod.residue is not None:
            mod_site_str = mod.residue
        mod_pos_str = mod.position if mod.position is not None else ''
        mod_site = ('%s%s' % (mod_site_str, mod_pos_str))
        site_states = states[mod.mod_type]
        if mod.is_modified:
            pattern[mod_site] = (site_states[1], WILD)
        else:
            pattern[mod_site] = (site_states[0], WILD)

    # Handle mutations
    for mc in agent.mutations:
        res_from = mc.residue_from if mc.residue_from else 'mut'
        res_to = mc.residue_to if mc.residue_to else 'X'
        if mc.position is None:
            mut_site_name = res_from
        else:
            mut_site_name = res_from + mc.position
        pattern[mut_site_name] = res_to

    # Handle location
    if agent.location is not None:
        pattern['loc'] = _n(agent.location)

    # Handle activity
    if agent.activity is not None:
        active_site_name = agent.activity.activity_type
        if agent.activity.is_active:
            active_site_state = 'active'
        else:
            active_site_state = 'inactive'
        pattern[active_site_name] = active_site_state

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

def set_extended_initial_condition(model, monomer=None, value=0):
    """Set an initial condition for monomers in "modified" state.

    This is useful when using downstream analysis that relies on reactions
    being active in the model. One example is BioNetGen-based reaction network
    diagram generation.
    """
    if not monomer:
        monomers = model.monomers
    else:
        monomers = [monomer]
    for monomer in monomers:
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
    url = get_identifiers_url(db_name, db_ref)
    if not url:
        return None
    subj = component
    ann = Annotation(subj, url, 'is')
    return ann

def parse_identifiers_url(url):
    """Parse an identifiers.org URL into (namespace, ID) tuple."""
    url_pattern = 'http://identifiers.org/([A-Za-z]+)/([A-Za-z0-9:]+)'
    match = re.match(url_pattern, url)
    if match is not None:
        g = match.groups()
        if not len(g) == 2:
            return (None, None)
        ns_map = {'hgnc': 'HGNC', 'uniprot': 'UP', 'chebi':'CHEBI',
                  'interpro':'IP', 'pfam':'XFAM'}
        ns = g[0]
        id = g[1]
        if not ns in ns_map.keys():
            return (None, None)
        if ns == 'hgnc':
            if id.startswith('HGNC:'):
                id = id[5:]
            else:
                logger.warning('HGNC URL missing "HGNC:" prefix: %s' % url)
                return (None, None)
        indra_ns = ns_map[ns]
        return (indra_ns, id)
    return (None, None)

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
        Example: one_step, interactions_only

        A dictionary of policies has keys corresponding to Statement types
        and values to the policy to be applied to that type of Statement.
        For Statement types whose policy is undefined, the 'default'
        policy is applied.
        Example: {'Phosphorylation': 'two_step'}

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
        self.default_initial_amount = 10000.0
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
        self.statements += stmts

    def make_model(self, policies=None, initial_conditions=True,
                   reverse_effects=False):
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
        ppa = PysbPreassembler(self.statements)
        ppa.replace_activities()
        if reverse_effects:
            ppa.add_reverse_effects()
        self.statements = ppa.statements
        # Set local policies for this make_model call that overwrite
        # the global policies of the PySB assembler
        if policies is not None:
            global_policies = self.policies
            if isinstance(policies, basestring):
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
            m.site_annotations = agent.site_annotations
            self.model.add_component(m)
            for db_name, db_ref in agent.db_refs.items():
                a = get_annotation(m, db_name, db_ref)
                if a is not None:
                    self.model.add_annotation(a)
            # Iterate over the active_forms
            for af in agent.active_forms:
                self.model.add_annotation(Annotation(m, af,
                                                     'has_active_pattern'))
            for iaf in agent.inactive_forms:
                self.model.add_annotation(Annotation(m, iaf,
                                                     'has_inactive_pattern'))
            for at in agent.activity_types:
                act_site_cond = {at: 'active'}
                self.model.add_annotation(Annotation(m, act_site_cond,
                                                     'has_active_pattern'))
                inact_site_cond = {at: 'inactive'}
                self.model.add_annotation(Annotation(m, inact_site_cond,
                                                     'has_inactive_pattern'))

        # Iterate over the statements to generate rules
        self._assemble()
        # Add initial conditions
        if initial_conditions:
            self.add_default_initial_conditions()

        # If local policies were applied, revert to the global one
        if policies is not None:
            self.policies = global_policies

        return self.model

    def add_default_initial_conditions(self, value=None):
        """Set default initial conditions in the PySB model.

        Parameters
        ----------
        value : Optional[float]
            Optionally a value can be supplied which will be the initial
            amount applied. Otherwise a built-in default is used.
        """
        if value is not None:
            try:
                value_num = float(value)
            except ValueError:
                logger.error('Invalid initial condition value.')
                return
        else:
            value_num = self.default_initial_amount
        if self.model is None:
            return
        for m in self.model.monomers:
            set_base_initial_condition(self.model, m, value_num)

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
        res = context_client.get_protein_expression(monomer_names, [cell_type])
        amounts = res.get(cell_type)
        if not amounts:
            logger.warning('Could not get context for %s cell type.' %
                           cell_type)
            self.add_default_initial_conditions()
            return
        monomers_found = []
        monomers_notfound = []
        for m in self.model.monomers:
            init = amounts.get(m.name)
            if init is not None:
                init_round = round(init)
                set_base_initial_condition(self.model, m, init_round)
                monomers_found.append(m.name)
            else:
                set_base_initial_condition(self.model, m,
                                           self.default_initial_amount)
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

    def export_model(self, format, file_name=None):
        """Save the assembled model in a modeling formalism other than PySB.

        For more details on exporting PySB models, see
        http://pysb.readthedocs.io/en/latest/modules/export/index.html

        Parameters
        ----------
        format : str
            The format to export into, for instance "kappa", "bngl",
            "sbml", "matlab", "mathematica", "potterswheel". See
            http://pysb.readthedocs.io/en/latest/modules/export/index.html
            for a list of supported formats. In addition to the formats
            supported by PySB itself, this method also provides "sbgn"
            output.

        file_name : Optional[str]
            An optional file name to save the exported model into.

        Returns
        -------
        exp_str : str
            The exported model string
        """
        # Handle SBGN as special case
        if format == 'sbgn':
            exp_str = export_sbgn(self.model)
        else:
            try:
                exp_str = pysb.export.export(self.model, format)
            except KeyError:
                logging.error('Unknown export format: %s' % format)
                return None

        if file_name:
            with open(file_name, 'wb') as fh:
                fh.write(exp_str.encode('utf-8'))
        return exp_str


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
            gene_mono.create_site(get_binding_site_name(bp))

complex_monomers_default = complex_monomers_one_step


def complex_assemble_one_step(stmt, model, agent_set):
    pairs = itertools.combinations(stmt.members, 2)
    for pair in pairs:
        agent1 = pair[0]
        agent2 = pair[1]
        param_name = agent1.name[0].lower() + \
                     agent2.name[0].lower() + '_bind'
        kf_bind = get_create_parameter(model, 'kf_' + param_name, 1e-6)
        kr_bind = get_create_parameter(model, 'kr_' + param_name, 1e-1)

        # Make a rule name
        rule_name = '_'.join([get_agent_rule_str(m) for m in pair])
        rule_name += '_bind'

        # Construct full patterns of each agent with conditions
        agent1_pattern = get_monomer_pattern(model, agent1)
        agent2_pattern = get_monomer_pattern(model, agent2)
        agent1_bs = get_binding_site_name(agent2)
        agent2_bs = get_binding_site_name(agent1)
        r = Rule(rule_name, agent1_pattern(**{agent1_bs: None}) + \
                            agent2_pattern(**{agent2_bs: None}) >>
                            agent1_pattern(**{agent1_bs: 1}) % \
                            agent2_pattern(**{agent2_bs: 1}),
                            kf_bind)
        anns = [Annotation(rule_name, agent1_pattern.monomer.name,
                           'rule_has_subject'),
                Annotation(rule_name, agent1_pattern.monomer.name,
                           'rule_has_object'),
                Annotation(rule_name, agent2_pattern.monomer.name,
                           'rule_has_subject'),
                Annotation(rule_name, agent2_pattern.monomer.name,
                           'rule_has_object'),
                Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
        add_rule_to_model(model, r, anns)


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
        anns = [Annotation(rule_name, monomer1_uncond.monomer.name,
                           'rule_has_subject'),
                Annotation(rule_name, monomer1_uncond.monomer.name,
                           'rule_has_object'),
                Annotation(rule_name, monomer2_uncond.monomer.name,
                           'rule_has_subject'),
                Annotation(rule_name, monomer2_uncond.monomer.name,
                           'rule_has_object'),
                Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
        add_rule_to_model(model, r, anns)


def complex_assemble_multi_way(stmt, model, agent_set):
    # Get the rate parameter
    abbr_name = ''.join([m.name[0].lower() for m in stmt.members])
    kf_bind = get_create_parameter(model, 'kf_' + abbr_name + '_bind', 1e-6)
    kr_bind = get_create_parameter(model, 'kr_' + abbr_name + '_bind', 1e-1)

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
            bp_bs = get_binding_site_name(bp)
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
            bound_bs = get_binding_site_name(bc.agent)
            gene_bs = get_binding_site_name(member)
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
    anns = [Annotation(rule_fwd.name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, rule_fwd, anns)
    anns = [Annotation(rule_rev.name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, rule_rev, anns)

complex_assemble_default = complex_assemble_one_step

# MODIFICATION ###################################################

def modification_monomers_interactions_only(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    act_type = mod_acttype_map[stmt.__class__]
    active_site = act_type
    enz.create_site(active_site)
    sub = agent_set.get_create_base_agent(stmt.sub)
    # See NOTE in monomers_one_step, below
    mod_condition_name = stmt.__class__.__name__.lower()
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))


def modification_monomers_one_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    # NOTE: This assumes that a Modification statement will only ever
    # involve a single phosphorylation site on the substrate (typically
    # if there is more than one site, they will be parsed into separate
    # Phosphorylation statements, i.e., phosphorylation is assumed to be
    # distributive. If this is not the case, this assumption will need to
    # be revisited.
    mod_condition_name = stmt.__class__.__name__.lower()
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))


def modification_monomers_two_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    mod_condition_name = stmt.__class__.__name__.lower()
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))

    # Create site for binding the substrate
    enz.create_site(get_binding_site_name(stmt.sub))
    sub.create_site(get_binding_site_name(stmt.enz))


def modification_assemble_interactions_only(stmt, model, agent_set):
    if stmt.enz is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    kr_bind = get_create_parameter(model, 'kr_bind', 1.0, unique=False)

    enz = model.monomers[stmt.enz.name]
    sub = model.monomers[stmt.sub.name]

    # See NOTE in monomers_one_step
    mod_condition_name = stmt.__class__.__name__.lower()
    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)

    rule_name = '%s_%s_%s_%s' % (rule_enz_str, mod_condition_name,
                                 rule_sub_str, mod_site)
    active_site = mod_acttype_map[stmt.__class__]
    # Create a rule specifying that the substrate binds to the kinase at
    # its active site
    lhs = enz(**{active_site: None}) + sub(**{mod_site: None})
    rhs = enz(**{active_site: 1}) % sub(**{mod_site: 1})
    r_fwd = Rule(rule_name + '_fwd', lhs >> rhs, kf_bind)
    add_rule_to_model(model, r_fwd)


def modification_assemble_one_step(stmt, model, agent_set, rate_law=None):
    if stmt.enz is None:
        return
    mod_condition_name = stmt.__class__.__name__.lower()

    # See NOTE in monomers_one_step
    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    rule_name = '%s_%s_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)

    # Remove pre-set activity flag
    enz_pattern = get_monomer_pattern(model, stmt.enz)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]
    sub_unmod = get_monomer_pattern(model, stmt.sub,
        extra_fields={mod_site: unmod_site_state})
    sub_mod = get_monomer_pattern(model, stmt.sub,
        extra_fields={mod_site: mod_site_state})


    if not rate_law:
        param_name = 'kf_%s%s_%s' % (stmt.enz.name[0].lower(),
                                      stmt.sub.name[0].lower(), mod_condition_name)
        mod_rate = get_create_parameter(model, param_name, 1e-6)
    elif rate_law == 'michaelis_menten':
        # Parameters
        param_name = ('Km_' + stmt.enz.name[0].lower() +
                      stmt.sub.name[0].lower() + '_' + mod_condition_name)
        Km = get_create_parameter(model, param_name, 1e8)
        param_name = ('kc_' + stmt.enz.name[0].lower() +
                      stmt.sub.name[0].lower() + '_' + mod_condition_name)
        kcat = get_create_parameter(model, param_name, 100)

        # We need an observable for the substrate to use in the rate law
        sub_obs = Observable(rule_name + '_sub_obs', sub_unmod)
        model.add_component(sub_obs)
        # Note that [E0]*[S] is automatically multiplied into this rate
        # as the product of the reactants therefore they don't appear
        # in this expression
        # v = Vmax*[S]/(Km+[S]) = kcat*[E0]*[S]/(Km + [S])
        mod_rate = Expression(rule_name + '_rate', kcat / (Km + sub_obs))
        model.add_component(mod_rate)

    r = Rule(rule_name,
            enz_pattern + sub_unmod >>
            enz_pattern + sub_mod,
            mod_rate)
    anns = [Annotation(rule_name, enz_pattern.monomer.name, 'rule_has_subject'),
            Annotation(rule_name, sub_unmod.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)


def modification_assemble_two_step(stmt, model, agent_set):
    mod_condition_name = stmt.__class__.__name__.lower()
    if stmt.enz is None:
        return
    sub_bs = get_binding_site_name(stmt.sub)
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
    kr_bind = get_create_parameter(model, param_name, 1e-1)
    param_name = ('kc_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_' + mod_condition_name)
    kf_mod = get_create_parameter(model, param_name, 100)

    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)

    enz_bs = get_binding_site_name(stmt.enz)
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]

    rule_name = '%s_%s_bind_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)
    r = Rule(rule_name,
        enz_unbound() + \
        sub_pattern(**{mod_site: unmod_site_state, enz_bs: None}) >>
        enz_bound() % \
        sub_pattern(**{mod_site: unmod_site_state, enz_bs: 1}),
        kf_bind)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

    rule_name = '%s_%s_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)
    r = Rule(rule_name,
        enz_bound() % \
            sub_pattern(**{mod_site: unmod_site_state, enz_bs: 1}) >>
        enz_unbound() + \
            sub_pattern(**{mod_site: mod_site_state, enz_bs: None}),
        kf_mod)
    anns = [Annotation(rule_name, enz_bound.monomer.name,
                       'rule_has_subject'),
            Annotation(rule_name, sub_pattern.monomer.name,
                       'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

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
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

modification_monomers_default = modification_monomers_one_step
modification_assemble_default = modification_assemble_one_step


# PHOSPHORYLATION ###################################################

def phosphorylation_monomers_atp_dependent(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    sub.create_mod_site(ist.ModCondition('phosphorylation',
                                         stmt.residue, stmt.position))
    # Create site for binding the substrate
    enz.create_site(get_binding_site_name(stmt.sub))
    sub.create_site(get_binding_site_name(stmt.enz))

    # Make ATP base agent and create binding sites
    atp = agent_set.get_create_base_agent(ist.Agent('ATP'))
    atp.create_site('b')
    enz.create_site('ATP')


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
    sub_bs = get_binding_site_name(stmt.sub)
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
    enz_bs = get_binding_site_name(stmt.enz)
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
    kr_bind_atp = get_create_parameter(model, param_name, 1.)
    rule_name = '%s_phospho_bind_atp' % (enz_rule_str)
    r = Rule(rule_name,
        enz_atp_unbound() + atp(b=None) >>
        enz_atp_bound() %  atp(b=1), kf_bind_atp)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

    # Enzyme releasing ATP
    rule_name = '%s_phospho_dissoc_atp' % (enz_rule_str)
    r = Rule(rule_name,
        enz_mon_uncond({atp_bs: 1}) % atp(b=1) >>
        enz_mon_uncond({atp_bs: None}) + atp(b=None), kr_bind_atp)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

    # Enzyme binding substrate
    param_name = ('kf_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kf_bind = get_create_parameter(model, param_name, 1e-6)
    param_name = ('kr_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_bind')
    kr_bind = get_create_parameter(model, param_name, 1e-1)
    param_name = ('kc_' + stmt.enz.name[0].lower() +
                  stmt.sub.name[0].lower() + '_phos')
    kf_phospho = get_create_parameter(model, param_name, 100)

    phos_site = get_mod_site_name('phosphorylation',
                                  stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    rule_name = '%s_phospho_bind_%s_%s' % \
        (rule_enz_str, rule_sub_str, phos_site)
    r = Rule(rule_name,
        enz_sub_unbound() + \
        sub_pattern(**{phos_site: 'u', enz_bs: None}) >>
        enz_sub_bound() % \
        sub_pattern(**{phos_site: 'u', enz_bs: 1}),
        kf_bind)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

    # Enzyme phosphorylating substrate
    rule_name = '%s_phospho_%s_%s' % \
        (rule_enz_str, rule_sub_str, phos_site)
    r = Rule(rule_name,
        enz_sub_atp_bound() % atp(b=2) % \
            sub_pattern(**{phos_site: 'u', enz_bs: 1}) >>
        enz_sub_atp_unbound() + atp(b=None) + \
            sub_pattern(**{phos_site: 'p', enz_bs: None}),
        kf_phospho)
    anns = [Annotation(rule_name, enz_sub_atp_bound.monomer.name,
                       'rule_has_subject'),
            Annotation(rule_name, sub_pattern.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)
    # Add rule annotations to model

    # Enzyme dissociating from substrate
    rule_name = '%s_dissoc_%s' % (enz_rule_str, sub_rule_str)
    r = Rule(rule_name, enz_mon_uncond(**{sub_bs: 1}) % \
             sub_mon_uncond(**{enz_bs: 1}) >>
             enz_mon_uncond(**{sub_bs: None}) + \
             sub_mon_uncond(**{enz_bs: None}), kr_bind)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)


# DEMODIFICATION #####################################################

def demodification_monomers_interactions_only(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    active_site = mod_acttype_map[stmt.__class__]
    enz.create_site(active_site)
    mod_condition_name = stmt.__class__.__name__.lower()[2:]
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))


def demodification_monomers_one_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    mod_condition_name = stmt.__class__.__name__.lower()[2:]
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))


def demodification_monomers_two_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    mod_condition_name = stmt.__class__.__name__.lower()[2:]
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))
    # Create site for binding the substrate
    enz.create_site(get_binding_site_name(stmt.sub))
    sub.create_site(get_binding_site_name(stmt.enz))


def demodification_assemble_interactions_only(stmt, model, agent_set):
    if stmt.enz is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    enz = model.monomers[stmt.enz.name]
    sub = model.monomers[stmt.sub.name]
    active_site = mod_acttype_map[stmt.__class__]
    # See NOTE in Phosphorylation.monomers_one_step
    demod_condition_name = stmt.__class__.__name__.lower()
    mod_condition_name = demod_condition_name[2:]
    demod_site = get_mod_site_name(mod_condition_name,
                                   stmt.residue, stmt.position)

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    r = Rule('%s_%s_%s_%s' %
             (rule_enz_str, demod_condition_name, rule_sub_str, demod_site),
             enz(**{active_site: None}) + sub(**{demod_site: None}) >>
             enz(**{active_site: 1}) % sub(**{demod_site: 1}),
             kf_bind)
    add_rule_to_model(model, r)


def demodification_assemble_one_step(stmt, model, agent_set, rate_law=None):
    if stmt.enz is None:
        return
    demod_condition_name = stmt.__class__.__name__.lower()
    mod_condition_name = demod_condition_name[2:]

    demod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    enz_pattern = get_monomer_pattern(model, stmt.enz)

    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]
    sub_unmod = get_monomer_pattern(model, stmt.sub,
        extra_fields={demod_site: unmod_site_state})
    sub_mod = get_monomer_pattern(model, stmt.sub,
        extra_fields={demod_site: mod_site_state})

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    rule_name = '%s_%s_%s_%s' % \
                (rule_enz_str, demod_condition_name, rule_sub_str, demod_site)

    if not rate_law:
        param_name = 'kf_%s%s_%s' % (stmt.enz.name[0].lower(),
                                     stmt.sub.name[0].lower(),
                                     demod_condition_name)
        mod_rate = get_create_parameter(model, param_name, 1e-6)
    elif rate_law == 'michaelis_menten':
        # Parameters
        param_name = ('Km_' + stmt.enz.name[0].lower() +
                      stmt.sub.name[0].lower() + '_' + mod_condition_name)
        Km = get_create_parameter(model, param_name, 1e8)
        param_name = ('kc_' + stmt.enz.name[0].lower() +
                      stmt.sub.name[0].lower() + '_' + mod_condition_name)
        kcat = get_create_parameter(model, param_name, 100)

        # We need an observable for the substrate to use in the rate law
        sub_obs = Observable(rule_name + '_sub_obs', sub_mod)
        model.add_component(sub_obs)
        mod_rate = Expression(rule_name + '_rate', kcat / (Km + sub_obs))
        model.add_component(mod_rate)

    r = Rule(rule_name,
             enz_pattern() + sub_mod >> enz_pattern() + sub_unmod,
             mod_rate)
    anns = [Annotation(r.name, enz_pattern.monomer.name, 'rule_has_subject'),
            Annotation(r.name, sub_mod.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)


def demodification_assemble_two_step(stmt, model, agent_set):
    if stmt.enz is None:
        return
    demod_condition_name = stmt.__class__.__name__.lower()
    mod_condition_name = demod_condition_name[2:]
    sub_bs = get_binding_site_name(stmt.sub)
    enz_bs = get_binding_site_name(stmt.enz)
    enz_bound = get_monomer_pattern(model, stmt.enz,
                                    extra_fields={sub_bs: 1})
    enz_unbound = get_monomer_pattern(model, stmt.enz,
                                      extra_fields={sub_bs: None})
    sub_pattern = get_monomer_pattern(model, stmt.sub)

    param_name = 'kf_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_bind'
    kf_bind = get_create_parameter(model, param_name, 1e-6)
    param_name = 'kr_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_bind'
    kr_bind = get_create_parameter(model, param_name, 1e-1)
    param_name = 'kc_' + stmt.enz.name[0].lower() + \
        stmt.sub.name[0].lower() + '_' + demod_condition_name
    kf_demod = get_create_parameter(model, param_name, 100)

    demod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    rule_name = '%s_%s_bind_%s_%s' % \
        (rule_enz_str, demod_condition_name, rule_sub_str, demod_site)
    r = Rule(rule_name,
             enz_unbound() + \
             sub_pattern(**{demod_site: mod_site_state, enz_bs: None}) >>
             enz_bound() % \
             sub_pattern(**{demod_site: mod_site_state, enz_bs: 1}),
             kf_bind)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

    rule_name = '%s_%s_%s_%s' % \
        (rule_enz_str, demod_condition_name, rule_sub_str, demod_site)
    r = Rule(rule_name,
        enz_bound() % \
            sub_pattern(**{demod_site: mod_site_state, enz_bs: 1}) >>
        enz_unbound() + \
            sub_pattern(**{demod_site: unmod_site_state, enz_bs: None}),
        kf_demod)
    anns = [Annotation(r.name, enz_bound.monomer.name, 'rule_has_subject'),
            Annotation(r.name, sub_pattern.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

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
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

demodification_monomers_default = demodification_monomers_one_step
demodification_assemble_default = demodification_assemble_one_step

# Map specific modification monomer/assembly functions to the generic
# Modification assembly function
policies = ['interactions_only', 'one_step', 'two_step', 'default']

mod_classes = [cls for cls in ist.AddModification.__subclasses__()]
for mc, func_type, pol in itertools.product(mod_classes,
                                            ('monomers', 'assemble'),
                                            policies):
    code = '{mc}_{func_type}_{pol} = ' \
            'modification_{func_type}_{pol}'.format(
                    mc=ist.modclass_to_modtype[mc], func_type=func_type,
                    pol=pol)
    exec(code)

demod_classes = [cls for cls in ist.RemoveModification.__subclasses__()]
for mc, func_type, pol in itertools.product(demod_classes,
                                            ('monomers', 'assemble'),
                                            policies):
    code = '{mc}_{func_type}_{pol} = ' \
            'demodification_{func_type}_{pol}'.format(
                    mc=ist.modclass_to_modtype[mc], func_type=func_type,
                    pol=pol)
    exec(code)

rate_laws = ['michaelis_menten']
for mc, rate_law in itertools.product(mod_classes, rate_laws):
    code = '{mc}_monomers_{rate_law} = {mc}_monomers_one_step'.format(
                mc=ist.modclass_to_modtype[mc], rate_law=rate_law)
    exec(code)
    code = '{mc}_assemble_{rate_law} = ' \
            'lambda a, b, c: modification_assemble_' \
            'one_step(a, b, c, "{rate_law}")'.format(
                mc=ist.modclass_to_modtype[mc], rate_law=rate_law)
    exec(code)

for mc, rate_law in itertools.product(demod_classes, rate_laws):
    code = '{mc}_monomers_{rate_law} = {mc}_monomers_one_step'.format(
                mc=ist.modclass_to_modtype[mc], rate_law=rate_law)
    exec(code)
    code = '{mc}_assemble_{rate_law} = ' \
            'lambda a, b, c: demodification_assemble_' \
            'one_step(a, b, c, "{rate_law}")'.format(
                    mc=ist.modclass_to_modtype[mc], rate_law=rate_law)
    exec(code)


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
    # http://www.jbc.org/content/286/4/2689.full
    kf_autophospho = get_create_parameter(model, param_name, 1e-2)

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
    anns = [Annotation(rule_name, pattern_unphos.monomer.name, 'rule_has_subject'),
            Annotation(rule_name, pattern_phos.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

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
    kf = get_create_parameter(model, param_name, 1e-6)

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
    anns = [Annotation(rule_name, enz_pattern.monomer.name, 'rule_has_subject'),
            Annotation(rule_name, sub_unphos.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

transphosphorylation_assemble_default = transphosphorylation_assemble_one_step

# ACTIVATION ######################################################

def regulateactivity_monomers_interactions_only(stmt, agent_set):
    subj = agent_set.get_create_base_agent(stmt.subj)
    obj = agent_set.get_create_base_agent(stmt.obj)
    if stmt.subj.activity is not None:
        subj_activity = stmt.subj.activity.activity_type
    else:
        subj_activity = 'activity'
    subj.create_site(subj_activity)
    obj.create_site(stmt.obj_activity)
    obj.create_site(stmt.obj_activity)


def regulateactivity_monomers_one_step(stmt, agent_set):
    subj = agent_set.get_create_base_agent(stmt.subj)
    obj = agent_set.get_create_base_agent(stmt.obj)
    # if stmt.subj_activity is not None:
    #    # Add the new active state flag to the list of active forms
    #    subj.add_activity_form({stmt.subj_activity: 'active'}, True)
    #    subj.add_activity_form({stmt.subj_activity: 'inactive'}, False)
    obj.create_site(stmt.obj_activity, ('inactive', 'active'))
    # Add the new active state flag to the list of active forms
    obj.add_activity_type(stmt.obj_activity)


def regulateactivity_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    subj = model.monomers[stmt.subj.name]
    obj = model.monomers[stmt.obj.name]

    if stmt.subj.activity:
        subj_activity = stmt.subj.activity.activity_type
    else:
        subj_activity = 'activity'

    subj_active_site = subj_activity
    obj_mod_site = stmt.obj_activity

    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_subj_str = get_agent_rule_str(stmt.subj)
    polarity_str = 'activates' if stmt.is_activation else 'deactivates'
    rule_name = '%s_%s_%s_%s' %\
             (rule_subj_str, polarity_str, rule_obj_str,
              stmt.obj_activity)
    r = Rule(rule_name,
             subj(**{subj_active_site: None}) +
             obj(**{obj_mod_site: None}) >>
             subj(**{subj_active_site: 1}) %
             obj(**{obj_mod_site: 1}),
             kf_bind)
    add_rule_to_model(model, r)


def regulateactivity_assemble_one_step(stmt, model, agent_set, rate_law=None):
    # This is the pattern coming directly from the subject Agent state
    # TODO: handle context here in conjunction with active forms
    subj_pattern = get_monomer_pattern(model, stmt.subj)

    obj_inactive = get_monomer_pattern(model, stmt.obj,
        extra_fields={stmt.obj_activity: 'inactive'})
    obj_active = get_monomer_pattern(model, stmt.obj,
        extra_fields={stmt.obj_activity: 'active'})

    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_subj_str = get_agent_rule_str(stmt.subj)
    polarity_str = 'activates' if stmt.is_activation else 'deactivates'
    rule_name = '%s_%s_%s_%s' % \
        (rule_subj_str, polarity_str, rule_obj_str,
         stmt.obj_activity)

    if not rate_law:
        param_name = 'kf_' + stmt.subj.name[0].lower() + \
            stmt.obj.name[0].lower() + '_act'
        act_rate = get_create_parameter(model, param_name, 1e-6)
    elif rate_law == 'michaelis_menten':
        # Parameters
        param_name = ('Km_' + stmt.subj.name[0].lower() +
                      stmt.obj.name[0].lower() + '_act')
        Km = get_create_parameter(model, param_name, 1e8)
        param_name = ('kc_' + stmt.subj.name[0].lower() +
                      stmt.obj.name[0].lower() + '_act')
        kcat = get_create_parameter(model, param_name, 100)

        # We need an observable for the substrate to use in the rate law
        obj_to_observe = obj_active if stmt.is_activation else obj_inactive
        obj_obs = Observable(rule_name + '_obj_obs', obj_to_observe)
        model.add_component(obj_obs)
        act_rate = Expression(rule_name + '_rate', kcat / (Km * obj_obs))
        model.add_component(act_rate) 

    obj_from, obj_to = (obj_inactive, obj_active) if stmt.is_activation else \
                       (obj_active, obj_inactive)

    r = Rule(rule_name, subj_pattern() + obj_from >> subj_pattern() + obj_to,
             act_rate)

    anns = [Annotation(rule_name, subj_pattern.monomer.name,
                       'rule_has_subject'),
            Annotation(rule_name, obj_active.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

regulateactivity_monomers_default = regulateactivity_monomers_one_step
regulateactivity_assemble_default = regulateactivity_assemble_one_step

activation_monomers_interactions_only = \
                    regulateactivity_monomers_interactions_only
activation_assemble_interactions_only = \
                    regulateactivity_assemble_interactions_only
activation_monomers_one_step = regulateactivity_monomers_one_step
activation_assemble_one_step = regulateactivity_assemble_one_step
activation_monomers_default = regulateactivity_monomers_one_step
activation_assemble_default = regulateactivity_assemble_one_step

activation_monomers_michaelis_menten = regulateactivity_monomers_one_step
activation_assemble_michaelis_menten = lambda a, b, c: \
    regulateactivity_assemble_one_step(a, b, c, 'michaelis_menten')

inhibition_monomers_interactions_only = \
                    regulateactivity_monomers_interactions_only
inhibition_assemble_interactions_only = \
                    regulateactivity_assemble_interactions_only
inhibition_monomers_one_step = regulateactivity_monomers_one_step
inhibition_assemble_one_step = regulateactivity_assemble_one_step
inhibition_monomers_default = regulateactivity_monomers_one_step
inhibition_assemble_default = regulateactivity_assemble_one_step

inhibition_monomers_michaelis_menten = regulateactivity_monomers_one_step
inhibition_assemble_michaelis_menten = lambda a, b, c: \
    regulateactivity_assemble_one_step(a, b, c, 'michaelis_menten')

# GEF #####################################################

def gef_monomers_interactions_only(stmt, agent_set):
    gef = agent_set.get_create_base_agent(stmt.gef)
    gef.create_site('gef_site')
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('p_loop')


def gef_monomers_one_step(stmt, agent_set):
    # Gef
    gef = agent_set.get_create_base_agent(stmt.gef)
    # Ras
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('gtpbound', ('inactive', 'active'))
    ras.add_activity_form({'gtpbound': 'active'}, True)
    ras.add_activity_form({'gtpbound': 'inactive'}, False)

gef_monomers_default = gef_monomers_one_step


def gef_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    gef = model.monomers[stmt.gef.name]
    ras = model.monomers[stmt.ras.name]
    rule_gef_str = get_agent_rule_str(stmt.gef)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_activates_%s' %
             (rule_gef_str, rule_ras_str),
             gef(**{'gef_site': None}) +
             ras(**{'p_loop': None}) >>
             gef(**{'gef_site': 1}) %
             ras(**{'p_loop': 1}),
             kf_bind)
    add_rule_to_model(model, r)


def gef_assemble_one_step(stmt, model, agent_set):
    gef_pattern = get_monomer_pattern(model, stmt.gef)
    ras_inactive = get_monomer_pattern(model, stmt.ras,
        extra_fields={'gtpbound': 'inactive'})
    ras_active = get_monomer_pattern(model, stmt.ras,
        extra_fields={'gtpbound': 'active'})

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
    anns = [Annotation(r.name, gef_pattern.monomer.name,
                       'rule_has_subject'),
            Annotation(r.name, ras_inactive.monomer.name, 'rule_has_object')]
    anns += [Annotation(r.name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

gef_assemble_default = gef_assemble_one_step

# GAP ####################################################

def gap_monomers_interactions_only(stmt, agent_set):
    gap = agent_set.get_create_base_agent(stmt.gap)
    gap.create_site('gap_site')
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('gtp_site')


def gap_monomers_one_step(stmt, agent_set):
    # Gap
    gap = agent_set.get_create_base_agent(stmt.gap)
    # Ras
    ras = agent_set.get_create_base_agent(stmt.ras)
    ras.create_site('gtpbound', ('inactive', 'active'))
    ras.add_activity_form({'gtpbound': 'active'}, True)
    ras.add_activity_form({'gtpbound': 'inactive'}, False)

gap_monomers_default = gap_monomers_one_step


def gap_assemble_interactions_only(stmt, model, agent_set):
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    gap = model.monomers[stmt.gap.name]
    ras = model.monomers[stmt.ras.name]
    rule_gap_str = get_agent_rule_str(stmt.gap)
    rule_ras_str = get_agent_rule_str(stmt.ras)
    r = Rule('%s_inactivates_%s' %
             (rule_gap_str, rule_ras_str),
             gap(**{'gap_site': None}) +
             ras(**{'gtp_site': None}) >>
             gap(**{'gap_site': 1}) %
             ras(**{'gtp_site': 1}),
             kf_bind)
    add_rule_to_model(model, r)


def gap_assemble_one_step(stmt, model, agent_set):
    gap_pattern = get_monomer_pattern(model, stmt.gap)
    ras_inactive = get_monomer_pattern(model, stmt.ras,
        extra_fields={'gtpbound': 'inactive'})
    ras_active = get_monomer_pattern(model, stmt.ras,
        extra_fields={'gtpbound': 'active'})

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
    anns = [Annotation(r.name, gap_pattern.monomer.name,
                       'rule_has_subject'),
            Annotation(r.name, ras_inactive.monomer.name, 'rule_has_object')]
    anns += [Annotation(r.name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

gap_assemble_default = gap_assemble_one_step

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

# GTPACTIVATION ######################################

gtpactivation_monomers_default = activation_monomers_default

gtpactivation_assemble_default = activation_assemble_default


# TRANSLOCATION ###############################################
def translocation_monomers_default(stmt, agent_set):
    # Skip if either from or to locations are missing
    if stmt.to_location is None:
        return
    agent = agent_set.get_create_base_agent(stmt.agent)
    from_loc = stmt.from_location if stmt.from_location else 'cytoplasm'
    states = [_n(from_loc), _n(stmt.to_location)]

    agent.create_site('loc', states)

def translocation_assemble_default(stmt, model, agent_set):
    if stmt.to_location is None:
        return
    from_loc = stmt.from_location if stmt.from_location else 'cytoplasm'
    param_name = 'kf_%s_%s_%s' % (_n(stmt.agent.name).lower(),
                                  _n(from_loc), _n(stmt.to_location))
    kf_trans = get_create_parameter(model, param_name, 1.0, unique=True)
    monomer = model.monomers[_n(stmt.agent.name)]
    rule_agent_str = get_agent_rule_str(stmt.agent)
    rule_name = '%s_translocates_%s_to_%s' % (rule_agent_str,
                                              _n(from_loc),
                                              _n(stmt.to_location))
    agent_from = get_monomer_pattern(model, stmt.agent,
                                     extra_fields={'loc': _n(from_loc)})
    agent_to = get_monomer_pattern(model, stmt.agent,
                                   extra_fields={'loc':
                                                 _n(stmt.to_location)})
    r = Rule(rule_name, agent_from >> agent_to, kf_trans)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

# DEGRADATION ###############################################

def decreaseamount_monomers_interactions_only(stmt, agent_set):
    if stmt.subj is None:
        return
    subj = agent_set.get_create_base_agent(stmt.subj)
    obj = agent_set.get_create_base_agent(stmt.obj)
    subj.create_site(get_binding_site_name(stmt.obj))
    obj.create_site(get_binding_site_name(stmt.subj))

def decreaseamount_monomers_one_step(stmt, agent_set):
    obj = agent_set.get_create_base_agent(stmt.obj)
    if stmt.subj is not None:
        subj = agent_set.get_create_base_agent(stmt.subj)

def decreaseamount_assemble_interactions_only(stmt, model, agent_set):
    # No interaction when subj is None
    if stmt.subj is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    subj_base_agent = agent_set.get_create_base_agent(stmt.subj)
    obj_base_agent = agent_set.get_create_base_agent(stmt.obj)
    subj = model.monomers[subj_base_agent.name]
    obj = model.monomers[obj_base_agent.name]
    rule_subj_str = get_agent_rule_str(stmt.subj)
    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_name = '%s_degrades_%s' % (rule_subj_str, rule_obj_str)

    subj_site_name = get_binding_site_name(stmt.obj)
    obj_site_name = get_binding_site_name(stmt.subj)

    r = Rule(rule_name,
             subj(**{subj_site_name: None}) + obj(**{obj_site_name: None}) >>
             subj(**{subj_site_name: 1}) % obj(**{obj_site_name: 1}),
             kf_bind)

    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    anns += [Annotation(rule_name, obj.name, 'rule_has_object')]
    if stmt.subj:
        anns += [Annotation(rule_name, subj.name, 'rule_has_subject')]
    add_rule_to_model(model, r, anns)

def decreaseamount_assemble_one_step(stmt, model, agent_set):
    obj_pattern = get_monomer_pattern(model, stmt.obj)
    rule_obj_str = get_agent_rule_str(stmt.obj)

    if stmt.subj is None:
        # See U. Alon paper on proteome dynamics at 10.1126/science.1199784 
        param_name = 'kf_' + stmt.obj.name[0].lower() + '_deg'
        kf_one_step_degrade = get_create_parameter(model, param_name, 2e-5,
                                                   unique=True)
        rule_name = '%s_degraded' % rule_obj_str
        r = Rule(rule_name, obj_pattern >> None, kf_one_step_degrade)
    else:
        subj_pattern = get_monomer_pattern(model, stmt.subj)
        # See U. Alon paper on proteome dynamics at 10.1126/science.1199784 
        param_name = 'kf_' + stmt.subj.name[0].lower() + \
                            stmt.obj.name[0].lower() + '_deg'
        # Scale the average apparent decreaseamount rate by the default
        # protein initial condition
        kf_one_step_degrade = get_create_parameter(model, param_name, 2e-9)
        rule_subj_str = get_agent_rule_str(stmt.subj)
        rule_name = '%s_degrades_%s' % (rule_subj_str, rule_obj_str)
        r = Rule(rule_name,
            subj_pattern + obj_pattern >> subj_pattern + None,
            kf_one_step_degrade)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    anns += [Annotation(rule_name, obj_pattern.monomer.name, 'rule_has_object')]
    if stmt.subj:
        anns += [Annotation(rule_name, subj_pattern.monomer.name, 'rule_has_subject')]
    add_rule_to_model(model, r, anns)

decreaseamount_assemble_default = decreaseamount_assemble_one_step
decreaseamount_monomers_default = decreaseamount_monomers_one_step

# SYNTHESIS ###############################################

increaseamount_monomers_interactions_only = \
                            decreaseamount_monomers_interactions_only

increaseamount_monomers_one_step = decreaseamount_monomers_one_step

def increaseamount_assemble_interactions_only(stmt, model, agent_set):
    # No interaction when subj is None
    if stmt.subj is None:
        return
    kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
    subj_base_agent = agent_set.get_create_base_agent(stmt.subj)
    obj_base_agent = agent_set.get_create_base_agent(stmt.obj)
    subj = model.monomers[subj_base_agent.name]
    obj = model.monomers[obj_base_agent.name]
    rule_subj_str = get_agent_rule_str(stmt.subj)
    rule_obj_str = get_agent_rule_str(stmt.obj)
    rule_name = '%s_synthesizes_%s' % (rule_subj_str, rule_obj_str)

    subj_site_name = get_binding_site_name(stmt.obj)
    obj_site_name = get_binding_site_name(stmt.subj)

    r = Rule(rule_name,
            subj(**{subj_site_name: None}) + obj(**{obj_site_name: None}) >>
            subj(**{subj_site_name: 1}) % obj(**{obj_site_name: 1}),
            kf_bind)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    anns += [Annotation(rule_name, obj.name, 'rule_has_object')]
    if stmt.subj:
        anns += [Annotation(rule_name, subj.name, 'rule_has_subject')]
    add_rule_to_model(model, r, anns)

def increaseamount_assemble_one_step(stmt, model, agent_set, rate_law=None):
    if stmt.subj is not None and (stmt.subj.name == stmt.obj.name):
        logger.warning('%s transcribes itself, skipping' % stmt.obj.name)
        return
    # We get the monomer pattern just to get a valid monomer
    # otherwise the patter will be replaced
    obj_pattern = get_monomer_pattern(model, stmt.obj)
    obj_monomer = obj_pattern.monomer
    # The obj Monomer needs to be synthesized in its "base" state
    # but it needs a fully specified monomer pattern
    sites_dict = {}
    for site in obj_monomer.sites:
        if site in obj_monomer.site_states:
            sites_dict[site] = obj_monomer.site_states[site][0]
        else:
            sites_dict[site] = None
    obj_pattern = obj_monomer(**sites_dict)
    rule_obj_str = get_agent_rule_str(stmt.obj)

    if stmt.subj is None:
        rule_name = '%s_synthesized' % rule_obj_str
        param_name = 'kf_' + stmt.obj.name[0].lower() + '_synth'
        kf_one_step_synth = get_create_parameter(model, param_name, 2,
                                                 unique=True)
        r = Rule(rule_name, None >> obj_pattern, kf_one_step_synth)
    else:
        subj_pattern = get_monomer_pattern(model, stmt.subj)
        rule_subj_str = get_agent_rule_str(stmt.subj)
        rule_name = '%s_synthesizes_%s' % (rule_subj_str, rule_obj_str)
        if not rate_law:
            param_name = 'kf_' + stmt.subj.name[0].lower() + \
                                stmt.obj.name[0].lower() + '_synth'
            # Scale the average apparent increaseamount rate by the default
            # protein initial condition
            synth_rate = get_create_parameter(model, param_name, 2e-4)
        if rate_law == 'hill':
            # k * [subj]**n / (K_A**n + [subj]**n)
            param_name = 'kf_' + stmt.subj.name[0].lower() + \
                                stmt.obj.name[0].lower() + '_synth'
            kf = get_create_parameter(model, param_name, 4)
            param_name = 'Ka_' + stmt.subj.name[0].lower() + \
                                stmt.obj.name[0].lower() + '_synth'
            Ka = get_create_parameter(model, param_name, 1e4)
            param_name = 'n_' + stmt.subj.name[0].lower() + \
                                stmt.obj.name[0].lower() + '_synth'
            n_hill = get_create_parameter(model, param_name, 1)
            subj_obs = Observable(rule_name + '_subj_obs', subj_pattern)
            model.add_component(subj_obs)
            synth_rate = Expression(rule_name + '_rate',
                kf * (subj_obs ** (n_hill-1)) / (Ka**n_hill + subj_obs**n_hill))
            model.add_component(synth_rate)

        r = Rule(rule_name, subj_pattern + None >> subj_pattern + obj_pattern,
                 synth_rate)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    anns += [Annotation(rule_name, obj_pattern.monomer.name, 'rule_has_object')]
    if stmt.subj:
        anns += [Annotation(rule_name, subj_pattern.monomer.name, 'rule_has_subject')]
    add_rule_to_model(model, r, anns)

increaseamount_monomers_default = increaseamount_monomers_one_step
increaseamount_assemble_default = increaseamount_assemble_one_step
increaseamount_monomers_hill = increaseamount_monomers_one_step
increaseamount_assemble_hill = lambda a, b, c: \
        increaseamount_assemble_one_step(a, b, c, 'hill')


# INFLUENCE ###################################################

influence_monomers_one_step = increaseamount_monomers_one_step
def influence_assemble_one_step(stmt, *args):
    if stmt.overall_polarity() == -1:
        return decreaseamount_assemble_one_step(stmt, *args)
    else:
        return increaseamount_assemble_one_step(stmt, *args, rate_law='hill')
influence_monomers_default = influence_monomers_one_step
influence_assemble_default = influence_assemble_one_step


# CONVERSION ###################################################

def conversion_monomers_one_step(stmt, agent_set):
    # Skip statements with more than one from object due to complications
    # with rate law
    if len(stmt.obj_from) > 1:
        return
    if stmt.subj is not None:
        agent_set.get_create_base_agent(stmt.subj)
    for obj in stmt.obj_from:
        agent_set.get_create_base_agent(obj)
    for obj in stmt.obj_to:
        agent_set.get_create_base_agent(obj)


def conversion_assemble_one_step(stmt, model, agent_set):
    # Skip statements with more than one from object due to complications
    # with rate law
    if len(stmt.obj_from) != 1:
        return

    # Create pieces needed for from object
    obj_from = stmt.obj_from[0]
    obj_from_pattern = get_monomer_pattern(model, obj_from)
    obj_from_monomer = obj_from_pattern.monomer
    rule_obj_from_str = get_agent_rule_str(obj_from)

    obj_to_monomers = [get_monomer_pattern(model, o).monomer for
                       o in stmt.obj_to]

    # Create pieces needed for to object
    # The obj Monomer needs to be synthesized in its "base" state
    # but it needs a fully specified monomer pattern
    obj_to_patterns = []
    for obj_to_monomer in obj_to_monomers:
        sites_dict = {}
        for site in obj_to_monomer.sites:
            if site in obj_to_monomer.site_states:
                sites_dict[site] = obj_to_monomer.site_states[site][0]
            else:
                sites_dict[site] = None
        obj_to_pattern = obj_to_monomer(**sites_dict)
        obj_to_patterns.append(obj_to_pattern)

    obj_to_pattern = ReactionPattern(obj_to_patterns)
    rule_obj_to_str = '_'.join([get_agent_rule_str(o) for o in stmt.obj_to])

    if stmt.subj is None:
        rule_name = '%s_converted_to_%s' % (rule_obj_from_str, rule_obj_to_str)
        param_name = 'kf_%s%s_convert' % (obj_from.name[0].lower(),
                                          obj_to_monomers[0].name[0].lower())
        kf_one_step_convert = get_create_parameter(model, param_name, 2,
                                                   unique=True)
        r = Rule(rule_name, obj_from_pattern >> obj_to_pattern,
                 kf_one_step_convert)
    else:
        subj_pattern = get_monomer_pattern(model, stmt.subj)
        result_pattern = obj_to_pattern
        result_pattern.complex_patterns.insert(0, subj_pattern)
        rule_subj_str = get_agent_rule_str(stmt.subj)
        rule_name = '%s_catalyzes_%s_converted_to_%s' % \
            (rule_subj_str, rule_obj_from_str, rule_obj_to_str)
        param_name = 'kf_%s%s%s_convert' % \
            (stmt.subj.name[0].lower(), obj_from.name[0].lower(),
             obj_to_monomers[0].name[0].lower())
        # Scale the average apparent increaseamount rate by the default
        # protein initial condition
        kf_one_step_convert = get_create_parameter(model, param_name, 2e-4)

        r = Rule(rule_name, subj_pattern + obj_from_pattern >>
                            result_pattern,
                 kf_one_step_convert)
    anns = [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)

conversion_monomers_default = conversion_monomers_one_step
conversion_assemble_default = conversion_assemble_one_step

class PysbPreassembler(object):
    def __init__(self, stmts=None):
        if not stmts:
            stmts = []
        self.statements = stmts
        self.agent_set = _BaseAgentSet()

    def add_statements(self, stmts):
        self.statements = stmts

    def _gather_active_forms(self):
        for stmt in self.statements:
            if isinstance(stmt, ist.ActiveForm):
                base_agent = self.agent_set.get_create_base_agent(stmt.agent)
                # Handle the case where an activity flag is set
                agent_to_add = stmt.agent
                if stmt.agent.activity:
                    new_agent = fast_deepcopy(stmt.agent)
                    new_agent.activity = None
                    agent_to_add = new_agent
                base_agent.add_activity_form(agent_to_add, stmt.is_active)

    def replace_activities(self):
        # TODO: handle activity hierarchies
        # First collect all explicit active forms
        self._gather_active_forms()
        new_stmts = []
        # Iterate over all statements
        for stmt in self.statements:
            stmt_agents = stmt.agent_list()
            num_agents = len(stmt_agents)
            # Make a list with an empty list for each Agent so that later
            # we can build combinations of Agent forms
            agent_forms = [[] for a in stmt_agents]
            for i, agent in enumerate(stmt_agents):
                # This is the case where there is an activity flag on an
                # Agent which we will attempt to replace with an explicit
                # active form
                if agent is not None and agent.activity is not None:
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
        # TODO: generalize to other modification sites
        pos_mod_sites = {}
        neg_mod_sites = {}
        syntheses = []
        degradations = []
        for stmt in self.statements:
            if isinstance(stmt, ist.Phosphorylation):
                agent = stmt.sub.name
                try:
                    pos_mod_sites[agent].append((stmt.residue, stmt.position))
                except KeyError:
                    pos_mod_sites[agent] = [(stmt.residue, stmt.position)]
            elif isinstance(stmt, ist.Dephosphorylation):
                agent = stmt.sub.name
                try:
                    neg_mod_sites[agent].append((stmt.residue, stmt.position))
                except KeyError:
                    neg_mod_sites[agent] = [(stmt.residue, stmt.position)]
            elif isinstance(stmt, ist.IncreaseAmount):
                syntheses.append(stmt.obj.name)
            elif isinstance(stmt, ist.DecreaseAmount):
                degradations.append(stmt.obj.name)

        new_stmts = []
        for agent_name, pos_sites in pos_mod_sites.items():
            neg_sites = neg_mod_sites.get(agent_name, [])
            no_neg_site = set(pos_sites).difference(set(neg_sites))
            for residue, position in no_neg_site:
                st = ist.Dephosphorylation(ist.Agent('phosphatase'),
                                           ist.Agent(agent_name),
                                           residue, position)
                new_stmts.append(st)
        for agent_name in syntheses:
            if agent_name not in degradations:
                st = ist.DecreaseAmount(None, ist.Agent(agent_name))
                new_stmts.append(st)

        self.statements += new_stmts

    @staticmethod
    def _set_agent_context(from_agent, to_agent):
        # TODO: what can we do about semantic conflicts here like the same
        # bound condition with True/False is_bound appearing in the
        # two contexts?
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
        to_agent.bound_conditions = \
            add_no_duplicate(to_agent.bound_conditions,
                             from_agent.bound_conditions)
        to_agent.mods = add_no_duplicate(to_agent.mods, from_agent.mods)
        to_agent.mutations = add_no_duplicate(to_agent.mutations,
                                              from_agent.mutations)
        to_agent.location = from_agent.location
        to_agent.activity = from_agent.activity

def export_sbgn(model):
    """Return an SBGN model string corresponding to the PySB model.

    This function first calls generate_equations on the PySB model to obtain
    a reaction network (i.e. individual species, reactions). It then iterates
    over each reaction and and instantiates its reactants, products, and the
    process itself as SBGN glyphs and arcs.

    Parameters
    ----------
    model : pysb.core.Model
        A PySB model to be exported into SBGN

    Returns
    -------
    sbgn_str : str
        An SBGN model as string
    """
    import lxml.etree
    import lxml.builder
    from pysb.bng import generate_equations
    from indra.assemblers.sbgn_assembler import SBGNAssembler

    logger.info('Generating reaction network with BNG for SBGN export. ' +
                'This could take a long time.')
    generate_equations(model)

    sa = SBGNAssembler()

    glyphs = {}
    for idx, species in enumerate(model.species):
        glyph = sa._glyph_for_complex_pattern(species)
        if glyph is None:
            continue
        sa._map.append(glyph)
        glyphs[idx] = glyph
    for reaction in model.reactions:
        # Get all the reactions / products / controllers of the reaction
        reactants = set(reaction['reactants']) - set(reaction['products'])
        products = set(reaction['products']) - set(reaction['reactants'])
        controllers = set(reaction['reactants']) & set(reaction['products'])
        # Add glyph for reaction
        process_glyph = sa._process_glyph('process')
        # Connect reactants with arcs
        for r in reactants:
            glyph = glyphs.get(r)
            if glyph is None:
                glyph_id = sa._none_glyph()
            else:
                glyph_id = glyph.attrib['id']
            sa._arc('consumption', glyph_id, process_glyph)
        # Connect products with arcs
        for p in products:
            glyph = glyphs.get(p)
            if glyph is None:
                glyph_id = sa._none_glyph()
            else:
                glyph_id = glyph.attrib['id']
            sa._arc('production', process_glyph, glyph_id)
        # Connect controllers with arcs
        for c in controllers:
            glyph = glyphs[c]
            sa._arc('catalysis', glyph.attrib['id'], process_glyph)

    sbgn_str = sa.print_model().decode('utf-8')
    return sbgn_str
