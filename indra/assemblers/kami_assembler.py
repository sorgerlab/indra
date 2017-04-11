from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import itertools
from collections import defaultdict
import indra.statements as ist
from indra.assemblers.pysb_assembler import \
        PysbAssembler, _is_whitelisted, \
        UnknownPolicyError, \
        get_binding_site_name, PysbPreassembler, \
        get_agent_rule_str, abbrevs, states, get_mod_site_name

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('kami_assembler')

class KamiAssembler(PysbAssembler):
    def make_model(self, policies=None, initial_conditions=True,
                   reverse_effects=False):
        """Assemble the Kami model from the collected INDRA Statements.

        This method assembles a Kami model from the set of INDRA Statements.
        The assembled model is both returned and set as the assembler's
        model argument.

        Parameters
        ----------
        policies : Optional[Union[str, dict]]
            A string or dictionary of policies, as defined in
            :py:class:`indra.assemblers.KamiAssembler`. This set of policies
            locally supersedes the default setting in the assembler. This
            is useful when this function is called multiple times with
            different policies.
        initial_conditions : Optional[bool]
            If True, default initial conditions are generated for the
            agents in the model.

        Returns
        -------
        model : dict
            The assembled Kami model.
        """
        ppa = PysbPreassembler(self.statements)
        ppa.replace_activities()
        if reverse_effects:
            ppa.add_reverse_effects()
        self.statements = ppa.statements
        # Set local policies for this make_model call that overwrite
        # the global policies of the Kami assembler
        if policies is not None:
            global_policies = self.policies
            if isinstance(policies, basestring):
                local_policies = {'other': policies}
            else:
                local_policies = {'other': 'default'}
                local_policies.update(policies)
            self.policies = local_policies

        self.model = {}
        graphs = []
        self.model['graphs'] = graphs
        self.model['typing'] = []

        # Action graph generated here
        action_graph = {'id': 'action_graph',
                        'attrs': {'name': 'action_graph'}}
        action_graph['graph'] = {'nodes': [], 'edges': []}
        graphs.append(action_graph)

        # Iterate over the statements to generate rules
        self._assemble()
        # Add initial conditions
        #if initial_conditions:
        #    self.add_default_initial_conditions()

        # If local policies were applied, revert to the global one
        if policies is not None:
            self.policies = global_policies

        return self.model

    def _assemble(self):
        """Calls the appropriate assemble method based on policies."""
        for stmt in self.statements:
            if _is_whitelisted(stmt):
                self._dispatch(stmt, 'assemble', self.model)

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
                #raise UnknownPolicyError('%s function %s not defined' %
                #                         (stage, func_name))
                logger.warning('%s function %s not defined' %
                               (stage, func_name))
                return
        return func(stmt, *args)


class Nugget(object):
    """Represents a Kami Nugget."""
    def __init__(self, id, name, rate):
        self.counters = defaultdict(int)
        self.id = id
        self.name = name
        self.rate = rate
        self.nodes = []
        self.edges = []
        self.typings = {}

    def add_agent(self, agent):
        """Add an INDRA Agent and its conditions to the Nugget."""
        agent_id = self.add_node(agent.name)
        self.add_typing(agent_id, 'agent')
        # Handle bound conditions
        for bc in agent.bound_conditions:
            # Here we make the assumption that the binding site
            # is simply named after the binding partner
            if bc.is_bound:
                test_type = 'is_bnd'
            else:
                test_type = 'is_free'
            bound_name = bc.agent.name
            agent_bs = get_binding_site_name(bc.agent)
            test_name = '%s_bound_to_%s_test' % (agent_id, bound_name)
            agent_bs_id = self.add_node(agent_bs)
            test_id = self.add_node(test_name)
            self.add_edge(agent_bs_id, agent_id)
            self.add_edge(agent_bs_id, test_id)
            self.add_typing(agent_bs_id, 'locus')
            self.add_typing(test_id, test_type)

        for mod in agent.mods:
            mod_site_str = abbrevs[mod.mod_type]
            if mod.residue is not None:
                mod_site_str = mod.residue
            mod_pos_str = mod.position if mod.position is not None else ''
            mod_site = ('%s%s' % (mod_site_str, mod_pos_str))
            site_states = states[mod.mod_type]
            if mod.is_modified:
                val = site_states[1]
            else:
                val = site_states[0]
            mod_site_id = self.add_node(mod_site, {'val': val})
            self.add_edge(mod_site_id, agent_id)
            self.add_typing(mod_site_id, 'state')
        return agent_id

    def add_node(self, name_base, attrs=None):
        """Add a node with a given base name to the Nugget and return ID."""
        if name_base not in self.counters:
            node_id = name_base
        else:
            node_id = '%s_%d' % (name_base, self.counters[name_base])
        node = {'id': node_id}
        if attrs:
            node['attrs'] = attrs
        self.nodes.append(node)
        self.counters[node_id] += 1
        return node_id

    def add_edge(self, from_node, to_node):
        """Add an edge between two nodes to the Nugget."""
        self.edges.append({'from': from_node, 'to': to_node})

    def add_typing(self, node_id, typing):
        """Add typing information to a node in the Nugget."""
        self.typings[node_id] = typing

    def get_nugget_dict(self):
        """Return the Nugget as a dictionary."""
        nugget_dict = \
            {'id': self.id,
             'graph': {
                 'nodes': self.nodes,
                 'edges': self.edges
                 },
             'attrs': {
                 'name': self.name,
                 'rate': self.rate
                 }
            }
        return nugget_dict

    def get_typing_dict(self):
        """Return the Nugget's typing information as a dictionary."""
        return self.typings


# COMPLEX ############################################################


def complex_assemble_one_step(stmt, model):
    pairs = itertools.combinations(stmt.members, 2)
    for pair in pairs:
        # Make a rule name
        nugget_name = '_'.join([get_agent_rule_str(m) for m in pair])
        nugget_name += '_bind'
        action_name =  nugget_name + '_act'
        kf_bind = 1e-6
        nugget = Nugget(nugget_name, nugget_name, kf_bind)
        action_id = nugget.add_node(action_name)
        # Initialize dicts/lists for this nugget
        nugget.add_typing(action_id, 'bnd')
        for agent in pair:
            agent_id = nugget.add_agent(agent)
            agent_bs = get_binding_site_name(agent)
            agent_bs_id = nugget.add_node(agent_bs)
            nugget.add_edge(agent_bs_id, agent_id)
            nugget.add_edge(agent_bs_id, action_id)
            # Add to the Kami typing dict
            nugget.add_typing(agent_bs_id, 'locus')
        # Typing dicts linking the nugget to the Action Graph and to the
        # Kami graph
        typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                          'mapping': {}, 'total': False,
                          'ignore_attrs': False}
        typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                            'mapping': nugget.get_typing_dict(), 'total': True,
                            'ignore_attrs': True}
        # Add the graphs for this nugget to the graphs and typing lists
        model['typing'] += [typing_dict_ag, typing_dict_kami]
        model['graphs'].append(nugget.get_nugget_dict())

        # In reverse reaction, assume that dissocition is unconditional
        nugget_name = '_'.join([get_agent_rule_str(m) for m in pair])
        nugget_name += '_dissociate'
        action_name =  nugget_name + '_act'
        kr_bind = 1e-1
        nugget = Nugget(nugget_name, nugget_name, kr_bind)
        action_id = nugget.add_node(action_name)
        nugget.add_typing(action_id, 'brk')
        for agent in pair:
            agent_bs = get_binding_site_name(agent)
            agent_id = nugget.add_node(agent.name)
            agent_bs_id = nugget.add_node(agent_bs)
            nugget.add_edge(agent_bs_id, agent_id)
            nugget.add_edge(agent_bs_id, action_id)
            nugget.add_typing(agent_id, 'agent')
            nugget.add_typing(agent_bs_id, 'locus')
        # Typing dicts linking the nugget to the Action Graph and to the
        # Kami graph
        typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                          'mapping': {}, 'total': False,
                          'ignore_attrs': False}
        typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                            'mapping': nugget.get_typing_dict(), 'total': True,
                            'ignore_attrs': True}
        # Add the graphs for this nugget to the graphs and typing lists
        model['typing'] += [typing_dict_ag, typing_dict_kami]
        model['graphs'].append(nugget.get_nugget_dict())


complex_assemble_default = complex_assemble_one_step


def _mod_demod_assemble_one_step(stmt, model, is_mod):
    # Define some basic parameters for the modification
    mod_condition_name = stmt.__class__.__name__.lower()
    if not is_mod:
        mod_condition_name = ist.modtype_to_inverse[mod_condition_name]

    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    nugget_name = '%s_%s_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)
    action_name =  nugget_name + '_act'
    kf_mod = 1e-6
    nugget = Nugget(nugget_name, nugget_name, kf_mod)
    enz_id = nugget.add_agent(stmt.enz)
    sub_id = nugget.add_agent(stmt.sub)

    st = states[mod_condition_name]
    from_state, to_state = (st[0], st[1]) if is_mod else (st[1], st[0])

    mod_site_id = nugget.add_node(mod_site, {'val': from_state})
    action_id = nugget.add_node(action_name, {'val': to_state})
    nugget.add_typing(mod_site_id, 'state')
    nugget.add_typing(action_id, 'mod')
    nugget.add_edge(mod_site_id, sub_id)
    nugget.add_edge(action_id, mod_site_id)
    nugget.add_edge(enz_id, action_id)

    # Typing dicts linking the nugget to the Action Graph and to the
    # Kami graph
    typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                      'mapping': {}, 'total': False,
                      'ignore_attrs': False}
    typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                        'mapping': nugget.get_typing_dict(), 'total': True,
                        'ignore_attrs': True}
    # Add the graphs for this nugget to the graphs and typing lists
    model['typing'] += [typing_dict_ag, typing_dict_kami]
    model['graphs'].append(nugget.get_nugget_dict())



def modification_assemble_one_step(stmt, model):
    if stmt.enz is None:
        return
    _mod_demod_assemble_one_step(stmt, model, True)


def demodification_assemble_one_step(stmt, model):
    if stmt.enz is None:
        return
    _mod_demod_assemble_one_step(stmt, model, False)


modification_assemble_default = modification_assemble_one_step
demodification_assemble_default = demodification_assemble_one_step


# Map specific modification monomer/assembly functions to the generic
# Modification assembly function
policies = ['one_step', 'default']

mod_classes = [cls for cls in ist.AddModification.__subclasses__()]
for mc, func_type, pol in itertools.product(mod_classes, ('assemble', ),
                                            policies):
    code = '{mc}_{func_type}_{pol} = ' \
            'modification_{func_type}_{pol}'.format(
                    mc=ist.modclass_to_modtype[mc], func_type=func_type,
                    pol=pol)
    exec(code)

demod_classes = [cls for cls in ist.RemoveModification.__subclasses__()]
for mc, func_type, pol in itertools.product(demod_classes, ('assemble', ),
                                            policies):
    code = '{mc}_{func_type}_{pol} = ' \
            'demodification_{func_type}_{pol}'.format(
                    mc=ist.modclass_to_modtype[mc], func_type=func_type,
                    pol=pol)
    exec(code)

