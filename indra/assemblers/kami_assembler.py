from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import itertools
import indra.statements as ist
from indra.assemblers.pysb_assembler import \
        PysbAssembler,\
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


# COMPLEX ############################################################

def complex_monomers_one_step(stmt, agent_set):
    pass


complex_monomers_default = complex_monomers_one_step


def complex_assemble_one_step(stmt, model, agent_set):
    pairs = itertools.combinations(stmt.members, 2)
    for pair in pairs:
        # Make a rule name
        nugget_name = '_'.join([get_agent_rule_str(m) for m in pair])
        nugget_name += '_bind'
        action_name =  nugget_name + '_act'
        kf_bind = 1e-6
        nugget_dict = {'id': nugget_name,
                       'graph': {},
                       'attrs': {'name': nugget_name,
                                  'rate': kf_bind}}
        # Initialize dicts/lists for this nugget
        nodes = [{'id': action_name}]
        edges = []
        typing_dict = {action_name: 'bnd'}
        for agent in pair:
            agent_nodes, agent_edges, agent_types = get_agent_conditions(agent)
            nodes += agent_nodes
            edges += agent_edges
            typing_dict.update(agent_types)
            # Construct full patterns of each agent with conditions
            #agent1_pattern = get_monomer_pattern(model, agent1)
            agent_bs = get_binding_site_name(agent)
            nodes.append({'id': agent.name})
            nodes.append({'id': agent_bs})
            edges.append({'from': agent_bs, 'to': agent.name})
            edges.append({'from': agent_bs, 'to': action_name})
            # Add to the Kami typing dict
            typing_dict.update({agent.name: 'agent', agent_bs: 'locus'})
        nugget_dict['graph']['nodes'] = nodes
        nugget_dict['graph']['edges'] = edges
        # Typing dicts linking the nugget to the Action Graph and to the
        # Kami graph
        typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                          'mapping': {}, 'total': False,
                          'ignore_attrs': False}
        typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                            'mapping': typing_dict, 'total': True,
                            'ignore_attrs': True}
        # Add the graphs for this nugget to the graphs and typing lists
        model['typing'] += [typing_dict_ag, typing_dict_kami]
        model['graphs'].append(nugget_dict)

        # In reverse reaction, assume that dissocition is unconditional
        nugget_name = '_'.join([get_agent_rule_str(m) for m in pair])
        nugget_name += '_dissociate'
        action_name =  nugget_name + '_act'
        kr_bind = 1e-1
        nugget_dict = {'id': nugget_name,
                       'attrs': {'name': nugget_name,
                                 'rate': kr_bind},
                       'graph': {}}
        # Initialize dicts/lists for this nugget
        nodes = [{'id': action_name}]
        edges = []
        typing_dict = {action_name: 'brk'}
        for agent in pair:
            agent_bs = get_binding_site_name(agent)
            nodes.append({'id': agent.name})
            nodes.append({'id': agent_bs})
            edges.append({'from': agent_bs, 'to': agent.name})
            edges.append({'from': agent_bs, 'to': action_name})
            # Add to the Kami typing dict
            typing_dict.update({agent.name: 'agent', agent_bs: 'locus'})
        nugget_dict['graph']['nodes'] = nodes
        nugget_dict['graph']['edges'] = edges
        # Typing dicts linking the nugget to the Action Graph and to the
        # Kami graph
        typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                          'mapping': {}, 'total': False,
                          'ignore_attrs': False}
        typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                            'mapping': typing_dict, 'total': True,
                            'ignore_attrs': True}
        # Add the graphs for this nugget to the graphs and typing lists
        model['typing'] += [typing_dict_ag, typing_dict_kami]
        model['graphs'].append(nugget_dict)


complex_assemble_default = complex_assemble_one_step



def modification_monomers_one_step(stmt, agent_set):
    pass

def modification_assemble_one_step(stmt, model, agent_set):
    if stmt.enz is None:
        return

    # Define some basic parameters for the modification
    mod_condition_name = stmt.__class__.__name__.lower()
    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]

    # Make a nugget name
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    nugget_name = '%s_%s_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)
    action_name =  nugget_name + '_act'
    kf_mod = 1e-6
    nugget_dict = {'id': nugget_name,
                   'graph': {},
                   'attrs': {'name': nugget_name,
                             'rate': kf_mod}}

    # Initialize dicts/lists for this nugget
    nodes = [{'id': action_name}]
    edges = []
    typing_dict = {action_name: 'bnd'}

    # Add enzyme conditions
    enz_nodes, enz_edges, enz_types = get_agent_conditions(stmt.enz)
    nodes += enz_nodes
    edges += enz_edges
    typing_dict.update(enz_types)

    # Add substrate conditions
    sub_nodes, sub_edges, sub_types = get_agent_conditions(stmt.sub)
    nodes += sub_nodes
    edges += sub_edges
    typing_dict.update(sub_types)

    # Add nodes/edges/types for the modification itself
    nodes.append({'id': mod_site, 'attrs': {'val': unmod_site_state}})
    nodes.append({'id': action_name, 'attrs': {'val': mod_site_state}})
    edges.append({'from': mod_site, 'to': stmt.sub.name})
    edges.append({'from': action_name, 'to': mod_site})
    edges.append({'from': stmt.enz.name, 'to': action_name})
    typing_dict.update({stmt.enz.name: 'agent', stmt.sub.name: 'agent',
                        mod_site: 'site', action_name: 'mod'})
    nugget_dict['graph']['nodes'] = nodes
    nugget_dict['graph']['edges'] = edges

    # Typing dicts linking the nugget to the Action Graph and to the
    # Kami graph
    typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                      'mapping': {}, 'total': False,
                      'ignore_attrs': False}
    typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                        'mapping': typing_dict, 'total': True,
                        'ignore_attrs': True}
    # Add the graphs for this nugget to the graphs and typing lists
    model['typing'] += [typing_dict_ag, typing_dict_kami]
    model['graphs'].append(nugget_dict)


def demodification_monomers_one_step(stmt, agent_set):
    pass


def demodification_assemble_one_step(stmt, model, agent_set):
    if stmt.enz is None:
        return

    # Define some basic parameters for the modification
    demod_condition_name = stmt.__class__.__name__.lower()
    mod_condition_name = ist.modtype_to_inverse[demod_condition_name]
    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]

    # Make a nugget name
    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    nugget_name = '%s_%s_%s_%s' % \
        (rule_enz_str, demod_condition_name, rule_sub_str, mod_site)
    action_name =  nugget_name + '_act'
    kf_mod = 1e-6
    nugget_dict = {'id': nugget_name,
                   'graph': {},
                   'attrs': {'name': nugget_name,
                             'rate': kf_mod}}

    # Initialize dicts/lists for this nugget
    nodes = [{'id': action_name}]
    edges = []
    typing_dict = {action_name: 'bnd'}

    # Add enzyme conditions
    enz_nodes, enz_edges, enz_types = get_agent_conditions(stmt.enz)
    nodes += enz_nodes
    edges += enz_edges
    typing_dict.update(enz_types)

    # Add substrate conditions
    sub_nodes, sub_edges, sub_types = get_agent_conditions(stmt.sub)
    nodes += sub_nodes
    edges += sub_edges
    typing_dict.update(sub_types)

    # Add nodes/edges/types for the modification itself
    nodes.append({'id': mod_site, 'attrs': {'val': mod_site_state}})
    nodes.append({'id': action_name, 'attrs': {'val': unmod_site_state}})
    edges.append({'from': mod_site, 'to': stmt.sub.name})
    edges.append({'from': action_name, 'to': mod_site})
    edges.append({'from': stmt.enz.name, 'to': action_name})
    typing_dict.update({stmt.enz.name: 'agent', stmt.sub.name: 'agent',
                        mod_site: 'site', action_name: 'mod'})
    nugget_dict['graph']['nodes'] = nodes
    nugget_dict['graph']['edges'] = edges

    # Typing dicts linking the nugget to the Action Graph and to the
    # Kami graph
    typing_dict_ag = {'from': nugget_name, 'to': 'action_graph',
                      'mapping': {}, 'total': False,
                      'ignore_attrs': False}
    typing_dict_kami = {'from': nugget_name, 'to': 'kami',
                        'mapping': typing_dict, 'total': True,
                        'ignore_attrs': True}
    # Add the graphs for this nugget to the graphs and typing lists
    model['typing'] += [typing_dict_ag, typing_dict_kami]
    model['graphs'].append(nugget_dict)


modification_monomers_default = modification_monomers_one_step
modification_assemble_default = modification_assemble_one_step
demodification_monomers_default = demodification_monomers_one_step
demodification_assemble_default = demodification_assemble_one_step


# Map specific modification monomer/assembly functions to the generic
# Modification assembly function
policies = ['one_step', 'default']

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


def get_agent_conditions(agent):
    nodes = []
    edges = []
    types = {}
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
        test_name = '%s_bound_to_%s_test' % (agent.name, bound_name)
        nodes += [{'id': agent_bs},
                  {'id': test_name}]
        edges += [{'from': agent_bs, 'to': agent.name},
                  {'from': agent_bs, 'to': test_name}]
        types.update({agent_bs: 'locus',
                      test_name: test_type})

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
        nodes += [{'id': mod_site, 'attrs': {'val': val}}]
        edges += [{'from': mod_site, 'to': agent.name}]
        types.update({mod_site: 'state'})

    return nodes, edges, types
    #TODO: locations, mutations, activity states
