from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.assemblers.pysb_assembler import \
        PysbAssembler,\
        UnknownPolicyError, \
        get_binding_site_name, PysbPreassembler, _BaseAgent, _BaseAgentSet, \
        get_agent_rule_str
import itertools

class KamiAssembler(PysbAssembler):
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
        self.model = {}
        graphs = []
        self.model['graphs'] = graphs
        self.model['typing'] = []
        self.agent_set = _BaseAgentSet()
        # Collect information about the monomers/self.agent_set from the
        # statements
        self._monomers()
        # Add the monomers to the model based on our BaseAgentSet
        # Action graph generated here
        action_graph = {'id': 'action_graph', 'name': 'action_graph'}
        action_graph['graph'] = {}
        graphs.append(action_graph)
        """
        for agent_name, agent in self.agent_set.items():
            nodes = []
            edges = []
            #m = Monomer(_n(agent_name), agent.sites, agent.site_states)
            #m.site_annotations = agent.site_annotations
            #self.model.add_component(m)
            #for db_name, db_ref in agent.db_refs.items():
            #    a = get_annotation(m, db_name, db_ref)
            #    if a is not None:
            #        self.model.add_annotation(a)
        """
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
        # Make a rule name
        rule_name = '_'.join([get_agent_rule_str(m) for m in pair])
        rule_name += '_bind'
        action_name =  rule_name + '_act'
        kf_bind = 1e-6
        kr_bind = 1e-1
        nugget_dict = {'id': rule_name, 'name': rule_name,
                       'graph': {'attributes':
                                    {'name': rule_name,
                                     'rate': kf_bind}}}
        # Initialize dicts/lists for this nugget
        nodes = [{'id': action_name}]
        edges = []
        typing_dict = {action_name: 'bnd'}
        for agent in pair:
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
        typing_dict_ag = {'from': rule_name, 'to': 'action_graph',
                          'typing': {}, 'total': False,
                          'ignore_attrs': False}
        typing_dict_kami = {'from': rule_name, 'to': 'kami',
                            'typing': typing_dict, 'total': True,
                            'ignore_attrs': True}
        # Add the graphs for this nugget to the graphs and typing lists
        model['typing'] += [typing_dict_ag, typing_dict_kami]
        model['graphs'].append(nugget_dict)

        """
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
        """

complex_assemble_default = complex_assemble_one_step




def phosphorylation_monomers_one_step(stmt, agent_set):
    if stmt.enz is None:
        return
    enz = agent_set.get_create_base_agent(stmt.enz)
    sub = agent_set.get_create_base_agent(stmt.sub)
    mod_condition_name = stmt.__class__.__name__.lower()
    sub.create_mod_site(ist.ModCondition(mod_condition_name,
                                         stmt.residue, stmt.position))

    """
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
    """

def phosphorylation_assemble_one_step(stmt, model, agent_set):
    """
    if stmt.enz is None:
        return
    mod_condition_name = stmt.__class__.__name__.lower()
    param_name = 'kf_%s%s_%s' % (stmt.enz.name[0].lower(),
                                  stmt.sub.name[0].lower(), mod_condition_name)
    kf_mod = get_create_parameter(model, param_name, 1e-6)

    # See NOTE in monomers_one_step
    mod_site = get_mod_site_name(mod_condition_name,
                                  stmt.residue, stmt.position)
    # Remove pre-set activity flag
    enz_pattern = get_monomer_pattern(model, stmt.enz)
    unmod_site_state = states[mod_condition_name][0]
    mod_site_state = states[mod_condition_name][1]
    sub_unmod = get_monomer_pattern(model, stmt.sub,
        extra_fields={mod_site: unmod_site_state})
    sub_mod = get_monomer_pattern(model, stmt.sub,
        extra_fields={mod_site: mod_site_state})

    rule_enz_str = get_agent_rule_str(stmt.enz)
    rule_sub_str = get_agent_rule_str(stmt.sub)
    rule_name = '%s_%s_%s_%s' % \
        (rule_enz_str, mod_condition_name, rule_sub_str, mod_site)
    r = Rule(rule_name,
            enz_pattern + sub_unmod >>
            enz_pattern + sub_mod,
            kf_mod)
    anns = [Annotation(rule_name, enz_pattern.monomer.name, 'rule_has_subject'),
            Annotation(rule_name, sub_unmod.monomer.name, 'rule_has_object')]
    anns += [Annotation(rule_name, stmt.uuid, 'from_indra_statement')]
    add_rule_to_model(model, r, anns)
    """
