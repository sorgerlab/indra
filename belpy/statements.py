import warnings
from sets import ImmutableSet
from pysb import *
from pysb import ReactionPattern, ComplexPattern, ComponentDuplicateNameError


abbrevs = {
    'PhosphorylationSerine': 'S',
    'PhosphorylationThreonine': 'T',
    'PhosphorylationTyrosine': 'Y',
    'Phosphorylation': 'phospho',
    'Ubiquitination': 'ub',
    'Farnesylation': 'farnesyl',
    'Hydroxylation': 'hydroxyl',
    'Acetylation': 'acetyl',
    'Sumoylation': 'sumo',
    'Glycosylation': 'glycosyl',
    'Methylation': 'methyl',
    'Modification': 'mod',
}

states = {
    'PhosphorylationSerine': ['u', 'p'],
    'PhosphorylationThreonine': ['u', 'p'],
    'PhosphorylationTyrosine': ['u', 'p'],
    'Phosphorylation': ['u', 'p'],
    'Ubiquitination': ['n', 'y'],
    'Farnesylation': ['n', 'y'],
    'Hydroxylation': ['n', 'y'],
    'Acetylation': ['n', 'y'],
    'Sumoylation': ['n', 'y'],
    'Glycosylation': ['n', 'y'],
    'Methylation': ['n', 'y'],
    'Modification': ['n', 'y'],
}

active_site_names = {
    'Kinase': 'kin_site',
    'Phosphatase': 'phos_site',
    'GtpBound': 'switch',
    'Catalytic': 'cat_site',
}

# The following dict specifies the default modification/binding site names for
# modifications resulting from a particular type of activity. For example, a
# protein with Kinase activity makes a modification of type "phospho" on its
# substrate, and a RasGTPase (with GtpBound activity) binds to a site of type
# "RBD" (Ras binding domain). This comes in handy for specifying
# ActivityActivity rules, where the modification site mediating the activation
# is not specified.
default_mod_site_names = {
    'Kinase': 'phospho',
    'GtpBound': 'RBD',
    'Phosphatase': 'phospho',
}

def site_name(stmt):
    """Return all site names for a modification-type statement."""
    names = []
    if isinstance(stmt.mod, (list, tuple)):
        for m, mp in zip(stmt.mod, stmt.mod_pos):
            mod = abbrevs[m]
            mod_pos = mp if mp is not None else ''
            names.append('%s%s' % (mod, mod_pos))
    else:
        mod = abbrevs[stmt.mod]
        mod_pos = stmt.mod_pos if stmt.mod_pos is not None else ''
        names.append('%s%s' % (mod, mod_pos))

    return names

def get_create_parameter(model, name, value):
    """Return parameter with given name, creating it if needed.

    If the parameter exists, the value is not changed; if it does not exist,
    the parameter value is set when it is created."""
    parameter = model.parameters.get(name)
    if parameter is None:
        parameter = Parameter(name, value)
        model.add_component(parameter)
    return parameter

class UnknownPolicyException(Exception):
    pass


class Agent(object):
    def __init__(self, name, mods=None, mod_sites=None, bound_to=None):
        self.name = name
        self.mods = mods
        self.mod_sites = mod_sites
        self.bound_to = bound_to


class Statement(object):
    """The parent class of all statements"""
    def __init__(self, stmt, citation, evidence, annotations):
        self.stmt = stmt
        self.citation = citation
        self.evidence = evidence
        self.annotations = annotations

    def monomers(self, agent_set, policies=None):
        """Calls the appropriate monomers method based on policies."""
        if policies is None or policies == 'one_step':
            self.monomers_one_step(agent_set)
        elif policies == 'interactions_only':
            self.monomers_interactions_only(agent_set)
        else:
            raise UnknownPolicyException(policies)

    def assemble(self, model, agent_set, policies=None):
        """Calls the appropriate assemble method based on policies."""
        if policies is None or policies == 'one_step':
            self.assemble_one_step(model, agent_set)
        elif policies == 'interactions_only':
            self.assemble_interactions_only(model, agent_set)
        else:
            raise UnknownPolicyException(policies)

    def monomers_one_step(self, agent_set):
        warnings.warn("%s.monomers_one_step not implemented" %
                      self.__class__.__name__)

    def assemble_one_step(self, model, agent_set):
        warnings.warn("%s.assemble_one_step not implemented" %
                      self.__class__.__name__)

    def monomers_interactions_only(self, agent_set):
        warnings.warn("%s.monomers_interactions_only not implemented" %
                      self.__class__.__name__)

    def assemble_interactions_only(self, model, agent_set):
        warnings.warn("%s.assemble_interactions_only not implemented" %
                      self.__class__.__name__)

class Modification(Statement):
    """Generic statement representing the modification of a protein"""
    def __init__(self, enz, sub, mod, mod_pos, stmt,
                 citation, evidence, annotations):
        super(Modification, self).__init__(stmt, citation, evidence,
                                           annotations)
        self.enz = enz
        self.sub = sub
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        return ("%s(%s, %s, %s, %s)" %
                (type(self).__name__, self.enz.name, self.sub.name, self.mod,
                 self.mod_pos))

class Phosphorylation(Modification):
    """Phosphorylation modification"""

    def monomers_interactions_only(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        enz.create_site(active_site_names['Kinase'])
        sub = agent_set.get_create_agent(self.sub.name)
        # See NOTE in monomers_one_step, below
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def monomers_one_step(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        enz.create_site('Kinase', ('inactive', 'active'))
        sub = agent_set.get_create_agent(self.sub.name)
        # NOTE: This assumes that a Phosphorylation statement will only ever
        # involve a single phosphorylation site on the substrate (typically
        # if there is more than one site, they will be parsed into separate
        # Phosphorylation statements, i.e., phosphorylation is assumed to be
        # distributive. If this is not the case, this assumption will need to
        # be revisited.
        sub.create_site(site_name(self)[0], ('u', 'p'))

        if self.enz.bound_to:
            enz_bound = agent_set.get_create_agent(self.enz.bound_to)
            enz_bound.create_site(self.enz.name)
            enz.create_site(self.enz_bound)

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.)

        enz = model.monomers[self.enz.name]
        sub = model.monomers[self.sub.name]

        # See NOTE in monomers_one_step
        site = site_name(self)[0]

        rule_name = '%s_phospho_%s_%s' % (self.enz.name, self.sub.name, site)
        active_site = active_site_names['Kinase']
        # Create a rule specifying that the substrate binds to the kinase at
        # its active site
        try:
            r = Rule(rule_name,
                     enz(**{active_site:None}) + sub(**{site:None}) <>
                     enz(**{active_site:1}) + sub(**{site:1}),
                     kf_bind, kf_bind)
            model.add_component(r)
        # If this rule is already in the model, issue a warning and continue
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

    def assemble_one_step(self, model, agent_set):
        kf_phospho = get_create_parameter(model, 'kf_phospho', 1e-6)

        enz = model.monomers[self.enz.name]
        sub = model.monomers[self.sub.name]

        # See NOTE in monomers_one_step
        site = site_name(self)[0]

        rule_name = '%s_phospho_%s_%s' % (self.enz.name, self.sub.name, site)
        try:
            # Iterate over all of the activating modification states of the
            # kinase
            enz_act_mods = agent_set[self.enz.name].activating_mods
            if enz_act_mods:
                for act_mod_pattern in enz_act_mods:
                    # Here we make the assumption that the binding site
                    # is simply named after the binding partner
                    if self.enz.bound_to:
                        act_mod_pattern[self.enz.bound_to] = 1
                        enz_bound = model.monomers[self.enz.bound_to]
                        enz_pattern = enz(**act_mod_pattern) % \
                                        enz_bound(**{self.enz.name:1})
                    else:
                        enz_pattern = enz(**act_mod_pattern)
                    r = Rule(rule_name,
                             enz_pattern + sub(**{site: 'u'}) >>
                             enz_pattern + sub(**{site: 'p'}),
                             kf_phospho)
                    model.add_component(r)
            # If there are no known activity modifications, we take this
            # statement as given and allow the enzyme to phosphorylate the
            # substrate unconditionally
            else:
                if self.enz.bound_to:
                    enz_bound = model.monomers[self.enz.bound_to]
                    enz_pattern = enz(**{self.enz.bound_to:1}) % \
                                    enz_bound(**{self.enz.name:1})
                else:
                    enz_pattern = enz()
                r = Rule(rule_name,
                         enz_pattern + sub(**{site: 'u'}) >>
                         enz_pattern + sub(**{site: 'p'}),
                         kf_phospho)
                model.add_component(r)

        # If this rule is already in the model, issue a warning and continue
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

class Hydroxylation(Modification):
    """Hydroxylation modification"""
    pass

class Sumoylation(Modification):
    """Sumoylation modification"""
    pass

class Acetylation(Modification):
    """Acetylation modification"""
    pass

class Ubiquitination(Modification):
    """Ubiquitination modification"""
    pass

class ActivityActivity(Statement):
    """Statement representing the activation of a protein as a result of the
    activity of another protein."""
    def __init__(self, subj, subj_activity, relationship, obj,
                 obj_activity, stmt, citation, evidence,
                 annotations):
        super(ActivityActivity, self).__init__(stmt,
                                               citation, evidence, annotations)
        self.subj = subj
        self.subj_activity = subj_activity
        self.obj = obj
        self.obj_activity = obj_activity
        self.relationship = relationship

    def monomers_interactions_only(self, agent_set):
        subj = agent_set.get_create_agent(self.subj.name)
        subj.create_site(active_site_names[self.subj_activity])
        obj = agent_set.get_create_agent(self.obj.name)
        obj.create_site(active_site_names[self.obj_activity])
        obj.create_site(default_mod_site_names[self.subj_activity])

    def assemble_interactions_only(self, model):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.)
        subj = model.monomers[self.subj.name]
        obj = model.monomers[self.obj.name]
        subj_active_site = active_site_names[self.subj_activity]
        obj_mod_site = default_mod_site_names[self.subj_activity]
        r = Rule('%s_%s_activates_%s_%s' %
                 (self.subj.name, self.subj_activity, self.obj.name,
                  self.obj_activity),
                 subj(**{subj_active_site: None}) +
                 obj(**{obj_mod_site: None}) >>
                 subj(**{subj_active_site: 1}) +
                 obj(**{obj_mod_site: 1}),
                 kf_bind)
        model.add_component(r)

    def monomers_one_step(self, model):
        subj = get_create_monomer(model, self.subj.name)
        create_site(subj, self.subj_activity, ('inactive', 'active'))
        obj = get_create_monomer(model, self.obj.name)
        create_site(obj, self.obj_activity, ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        kf_one_step_activate = \
                       get_create_parameter(model, 'kf_one_step_activate', 1e-6)

        subj = model.monomers[self.subj.name]
        obj = model.monomers[self.obj.name]

        r = Rule('%s_%s_activates_%s_%s' %
                 (self.subj.name, self.subj_activity, self.obj.name,
                  self.obj_activity),
                 subj(**{self.subj_activity: 'active'}) +
                 obj(**{self.obj_activity: 'inactive'}) >>
                 subj(**{self.subj_activity: 'active'}) +
                 obj(**{self.obj_activity: 'active'}),
                 kf_one_step_activate)
        model.add_component(r)

    def __str__(self):
        return ("%s(%s, %s, %s, %s, %s)" %
                (type(self).__name__, self.subj.name, self.subj_activity,
                 self.relationship, self.obj.name, self.obj_activity))

class RasGtpActivityActivity(ActivityActivity):
    pass

class Dephosphorylation(Statement):
    def __init__(self, phos, sub, mod, mod_pos, stmt,
                 citation, evidence, annotations):
        super(Dephosphorylation, self).__init__(stmt, citation,
                                                evidence, annotations)
        self.phos = phos
        self.sub = sub
        self.mod = mod
        self.mod_pos = mod_pos

    def monomers_interactions_only(self, agent_set):
        phos = agent_set.get_create_agent(self.phos.name)
        phos.create_site(active_site_names['Phosphatase'])
        sub = agent_set.get_create_agent(self.sub.name)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.)
        phos = model.monomers[self.phos.name]
        sub = model.monomers[self.sub.name]
        phos_site = active_site_names['Phosphatase']
        # See NOTE in Phosphorylation.monomers_one_step
        site = site_name(self)[0]
        r = Rule('%s_dephospho_%s_%s' %
                 (self.phos.name, self.sub.name, site),
                 phos(**{phos_site: None}) + sub(**{site: None}) >>
                 phos(**{phos_site: 1}) + sub(**{site: 1}),
                 kf_bind)
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        phos = agent_set.get_create_agent(self.phos.name)
        sub = agent_set.get_create_agent(self.sub.name)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_one_step(self, model, agent_set):
        kf_dephospho = get_create_parameter(model, 'kf_dephospho', 1e-6)

        phos = model.monomers[self.phos.name]
        sub = model.monomers[self.sub.name]

        site = site_name(self)[0]
        r = Rule('%s_dephospho_%s_%s' %
                 (self.phos.name, self.sub.name, site),
                 phos() + sub(**{site: 'p'}) >>
                 phos() + sub(**{site: 'u'}),
                 kf_dephospho)
        model.add_component(r)

    def __str__(self):
        return ("Dephosphorylation(%s, %s, %s, %s)" %
                (self.phos.name, self.sub.name, self.mod, self.mod_pos))

class ActivityModification(Statement):
    """Statement representing the activation of a protein as a result
    of a residue modification"""
    def __init__(self, monomer, mod, mod_pos, relationship, activity,
                 stmt, citation, evidence, annotations):
        super(ActivityModification, self).__init__(stmt, citation,
                                                   evidence, annotations)
        self.monomer = monomer
        self.mod = mod
        self.mod_pos = mod_pos
        self.relationship = relationship
        self.activity = activity

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def monomers_one_step(self, agent_set):
        agent = agent_set.get_create_agent(self.monomer.name)
        sites = site_name(self)
        active_states = [states[m][1] for m in self.mod]

        activity_pattern = {}
        for i, s in enumerate(sites):
            site_states = states[self.mod[i]]
            active_state = site_states[1]
            agent.create_site(s, site_states)
            activity_pattern[s] = active_state
        # Add the site/state for the activity itself FIXME FIXME FIXME
        agent.create_site(self.activity, ('inactive', 'active'))

        # Add this activity modification explicitly to the agent's list
        # of activating modifications
        agent.add_activating_modification(activity_pattern)
        # Inactivating modifications will require a different treatment
        # of the resolution of when the agent is active
        if self.relationship == 'DirectlyDecreases':
            warnings.warn('Inactivating modifications not currently '
                          'implemented!')

    def assemble_one_step(self, model, agent_set):
        pass

    def __str__(self):
        return ("ActivityModification(%s, %s, %s, %s, %s)" %
                (self.monomer.name, self.mod, self.mod_pos, self.relationship,
                 self.activity))

class ActivatingSubstitution(Statement):
    """Statement representing the activation of a protein as a result
    of a residue substitution"""
    def __init__(self, monomer, wt_residue, pos, sub_residue, activity,
                 stmt, citation, evidence, annotations):
        super(ActivatingSubstitution, self).__init__(stmt, citation,
                                                     evidence, annotations)
        self.monomer.name = monomer.name
        self.wt_residue = wt_residue
        self.pos = pos
        self.sub_residue = sub_residue
        self.activity = activity

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def __str__(self):
        return ("ActivatingSubstitution(%s, %s, %s, %s, %s)" %
                (self.monomer.name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity))

class RasGef(Statement):
    """Statement representing the activation of a GTP-bound protein
    upon Gef activity."""

    def __init__(self, gef, gef_activity, ras,
                 stmt, citation, evidence, annotations):
        super(RasGef, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gef = gef
        self.gef_activity = gef_activity
        self.ras = ras

    def monomers_interactions_only(self, agent_set):
        gef = agent_set.get_create_agent(self.gef.name)
        gef.create_site('gef_site')
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('p_loop')

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.)
        gef = model.monomers[self.gef.name]
        ras = model.monomers[self.ras.name]
        r = Rule('%s_activates_%s' %
                 (self.gef.name, self.ras.name),
                 gef(**{'gef_site':None}) +
                 ras(**{'p_loop':None}) >>
                 gef(**{'gef_site': 1}) +
                 ras(**{'p_loop': 1}),
                 kf_bind)
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        gef = agent_set.get_create_agent(self.gef.name)
        gef.create_site(self.gef_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        kf_gef = get_create_parameter(model, 'kf_gef', 1e-6)

        gef = model.monomers[self.gef.name]
        ras = model.monomers[self.ras.name]

        r = Rule('%s_activates_%s' %
                 (self.gef.name, self.ras.name),
                 gef(**{self.gef_activity: 'active'}) +
                 ras(**{'GtpBound': 'inactive'}) >>
                 gef(**{self.gef_activity: 'active'}) +
                 ras(**{'GtpBound': 'active'}),
                 kf_gef)
        model.add_component(r)

    def __str__(self):
        return ("RasGef(%s, %s, %s)" %
                (self.gef.name, self.gef_activity, self.ras.name))


class RasGap(Statement):
    """Statement representing the inactivation of a GTP-bound protein
    upon Gap activity."""
    def __init__(self, gap, gap_activity, ras,
                 stmt, citation, evidence, annotations):
        super(RasGap, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gap = gap
        self.gap_activity = gap_activity
        self.ras = ras

    def monomers_interactions_only(self, agent_set):
        gap = agent_set.get_create_agent(self.gap.name)
        gap.create_site('gap_site')
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('gtp_site')

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.)
        gap = model.monomers[self.gap.name]
        ras = model.monomers[self.ras.name]
        r = Rule('%s_inactivates_%s' %
                 (self.gap.name, self.ras.name),
                 gap(**{'gap_site': None}) +
                 ras(**{'gtp_site': None}) >>
                 gap(**{'gap_site': 1}) +
                 ras(**{'gtp_site': 1}),
                 kf_bind)
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        gap = agent_set.get_create_agent(self.gap.name)
        gap.create_site(self.gap_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        kf_gap = get_create_parameter(model, 'kf_gap', 1e-6)

        gap = model.monomers[self.gap.name]
        ras = model.monomers[self.ras.name]

        r = Rule('%s_inactivates_%s' %
                 (self.gap.name, self.ras.name),
                 gap(**{self.gap_activity: 'active'}) +
                 ras(**{'GtpBound': 'active'}) >>
                 gap(**{self.gap_activity: 'active'}) +
                 ras(**{'GtpBound': 'inactive'}),
                 kf_gap)
        model.add_component(r)

    def __str__(self):
        return ("RasGap(%s, %s, %s)" %
                (self.gap.name, self.gap_activity, self.ras.name))


class Complex(Statement):
    """Statement representing complex formation between a set of members"""
    def __init__(self, members):
        self.members = members

    def monomers_interactions_only(self, agent_set):
        return self.monomers_one_step(agent_set)

    def assemble_interactions_only(self, model, agent_set):
        return self.assemble_one_step(model, agent_set)

    def monomers_one_step(self, agent_set):
        """In this (very simple) implementation, proteins in a complex are
        each given site names corresponding to each of the other members
        of the complex. So the resulting complex is "fully connected" in
        that each is specified as bound to all the others."""
        for i, member in enumerate(self.members):
            gene_name = member.name
            gene_mono = agent_set.get_create_agent(gene_name)
            if member.bound_to:
                bound_name = member.bound_to
                bound_mono = agent_set.get_create_agent(bound_name)
                gene_mono.create_site(bound_name)
                bound_mono.create_site(gene_name)
            # Specify a binding site for each of the other complex members
            # bp = abbreviation for "binding partner"
            for j, bp in enumerate(self.members):
                # The protein doesn't bind to itself!
                if i == j:
                    continue
                gene_mono.create_site(bp.name)

    def assemble_one_step(self, model, agent_set):
        # Get the rate parameter
        kf_bind = get_create_parameter(model, 'kf_bind', 1e-6)

        # Make a rule name
        rule_name = '_'.join([m.name for m in self.members])
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
        for i, member in enumerate(self.members):
            gene_name = member.name
            mono = model.monomers[gene_name]
            # Specify free and bound states for binding sites for each of
            # the other complex members
            # (bp = abbreviation for "binding partner")
            left_site_dict = {}
            right_site_dict = {}
            for j, bp in enumerate(self.members):
                # The protein doesn't bind to itself!
                if i == j:
                    continue
                # Check to see if we've already created a bond index for these
                # two binding partners
                bp_set = ImmutableSet([i, j])
                if bp_set in bond_indices:
                    bond_ix = bond_indices[bp_set]
                # If we haven't see this pair of proteins yet, add a new bond
                # index to the dict
                else:
                    bond_ix = bond_counter
                    bond_indices[bp_set] = bond_ix
                    bond_counter += 1
                # Fill in the entries for the site dicts
                left_site_dict[bp.name] = None
                right_site_dict[bp.name] = bond_ix
            if member.bound_to:
                bound_name = member.bound_to
                bound = model.monomers[bound_name]
                left_site_dict[bound_name] = bond_counter
                right_site_dict[bound_name] = bond_counter
                left_pattern = mono(**left_site_dict) % \
                               bound(**{gene_name:bond_counter})
                right_pattern = mono(**right_site_dict) % \
                                bound(**{gene_name:bond_counter})
                bond_counter += 1
            else:
                left_pattern = mono(**left_site_dict)
                right_pattern = mono(**right_site_dict)
            # Build up the left- and right-hand sides of the rule from
            # monomer patterns with the appropriate site dicts
            lhs = lhs + left_pattern
            rhs = rhs % right_pattern
        # Finally, create the rule and add it to the model
        try:
            rule = Rule(rule_name, lhs <> rhs, kf_bind, kf_bind)
            model.add_component(rule)
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

    def __str__(self):
        return ("Complex(%s)" % [m.name for m in self.members])
