import warnings
from pysb import *
from pysb import ReactionPattern, ComplexPattern, ComponentDuplicateNameError
from sets import ImmutableSet

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

class Statement(object):
    """The parent class of all statements"""
    def __init__(self, stmt, citation, evidence, annotations):
        self.stmt = stmt
        self.citation = citation
        self.evidence = evidence
        self.annotations = annotations

    def monomers(self, model):
        warnings.warn("%s.monomers not implemented" % self.__class__.__name__)

    def assemble(self, model):
        warnings.warn("%s.assemble not implemented" % self.__class__.__name__)

class Modification(Statement):
    """Generic statement representing the modification of a protein"""
    def __init__(self, enz_name, sub_name, mod, mod_pos, stmt,
                 citation, evidence, annotations):
        super(Modification, self).__init__(stmt, citation, evidence,
                                           annotations)
        self.enz_name = enz_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        return ("%s(%s, %s, %s, %s)" %
                (type(self).__name__, self.enz_name, self.sub_name, self.mod,
                 self.mod_pos))

class Phosphorylation(Modification):
    """Phosphorylation modification"""
    def monomers(self, agent_set):
        enz = agent_set.get_create_agent(self.enz_name)
        enz.create_site('Kinase', ('inactive', 'active'))
        sub = agent_set.get_create_agent(self.sub_name)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble(self, model, agent_set):
        try:
            kf_phospho = model.parameters['kf_phospho']
        except KeyError:
            kf_phospho = Parameter('kf_phospho', 1e-6)
            model.add_component(kf_phospho)

        enz = model.monomers[self.enz_name]
        sub = model.monomers[self.sub_name]

        site = site_name(self)[0]

        rule_name = '%s_phospho_%s_%s' % (self.enz_name, self.sub_name, site)
        try:
            # Iterate over all of the activating modification states of the
            # kinase
            enz_act_mods = agent_set[self.enz_name].activating_mods
            if enz_act_mods:
                for act_mod_pattern in enz_act_mods:
                    r = Rule(rule_name,
                             enz(**act_mod_pattern) + sub(**{site: 'u'}) >>
                             enz(**act_mod_pattern) + sub(**{site: 'p'}),
                             kf_phospho)
                    model.add_component(r)
            # If there are no known activity modifications, we take this
            # statement as given and allow the enzyme to phosphorylate the
            # substrate unconditionally
            else:
                r = Rule(rule_name,
                         enz() + sub(**{site: 'u'}) >>
                         enz() + sub(**{site: 'p'}),
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
    def __init__(self, subj_name, subj_activity, relationship, obj_name,
                 obj_activity, stmt, citation, evidence,
                 annotations):
        super(ActivityActivity, self).__init__(stmt,
                                               citation, evidence, annotations)
        self.subj_name = subj_name
        self.subj_activity = subj_activity
        self.obj_name = obj_name
        self.obj_activity = obj_activity
        self.relationship = relationship

    def monomers(self, model):
        subj = get_create_monomer(model, self.subj_name)
        create_site(subj, self.subj_activity, ('inactive', 'active'))
        obj = get_create_monomer(model, self.obj_name)
        create_site(obj, self.obj_activity, ('inactive', 'active'))

    def assemble(self, model):
        try:
            kf_one_step_activate = model.parameters['kf_one_step_activate']
        except KeyError:
            kf_one_step_activate = Parameter('kf_one_step_activate', 1e-6)
            model.add_component(kf_one_step_activate)

        subj = model.monomers[self.subj_name]
        obj = model.monomers[self.obj_name]

        r = Rule('%s_%s_activates_%s_%s' %
                 (self.subj_name, self.subj_activity, self.obj_name,
                  self.obj_activity),
                 subj(**{self.subj_activity: 'active'}) +
                 obj(**{self.obj_activity: 'inactive'}) >>
                 subj(**{self.subj_activity: 'active'}) +
                 obj(**{self.obj_activity: 'active'}),
                 kf_one_step_activate)
        model.add_component(r)

    def __str__(self):
        return ("%s(%s, %s, %s, %s, %s)" %
                (type(self).__name__, self.subj_name, self.subj_activity,
                 self.relationship, self.obj_name, self.obj_activity))

class RasGtpActivityActivity(ActivityActivity):
    pass

class Dephosphorylation(Statement):
    def __init__(self, phos_name, sub_name, mod, mod_pos, stmt,
                 citation, evidence, annotations):
        super(Dephosphorylation, self).__init__(stmt, citation,
                                                evidence, annotations)
        self.phos_name = phos_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def monomers(self, agent_set):
        phos = agent_set.get_create_agent(self.phos_name)
        sub = agent_set.get_create_agent(self.sub_name)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble(self, model, agent_set):
        try:
            kf_dephospho = model.parameters['kf_dephospho']
        except KeyError:
            kf_dephospho = Parameter('kf_dephospho', 1e-6)
            model.add_component(kf_dephospho)

        phos = model.monomers[self.phos_name]
        sub = model.monomers[self.sub_name]

        site = site_name(self)[0]
        r = Rule('%s_dephospho_%s_%s' %
                 (self.phos_name, self.sub_name, site),
                 phos() + sub(**{site: 'p'}) >>
                 phos() + sub(**{site: 'u'}),
                 kf_dephospho)
        model.add_component(r)

    def __str__(self):
        return ("Dehosphorylation(%s, %s, %s, %s)" %
                (self.phos_name, self.sub_name, self.mod, self.mod_pos))

class ActivityModification(Statement):
    """Statement representing the activation of a protein as a result
    of a residue modification"""
    def __init__(self, monomer_name, mod, mod_pos, relationship, activity,
                 stmt, citation, evidence, annotations):
        super(ActivityModification, self).__init__(stmt, citation,
                                                   evidence, annotations)
        self.monomer_name = monomer_name
        self.mod = mod
        self.mod_pos = mod_pos
        self.relationship = relationship
        self.activity = activity

    def monomers(self, agent_set):
        agent = agent_set.get_create_agent(self.monomer_name)
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

    def assemble(self, model, agent_set):
        pass

    def __str__(self):
        return ("ActivityModification(%s, %s, %s, %s, %s)" %
                (self.monomer_name, self.mod, self.mod_pos, self.relationship,
                 self.activity))

class ActivatingSubstitution(Statement):
    """Statement representing the activation of a protein as a result
    of a residue substitution"""
    def __init__(self, monomer_name, wt_residue, pos, sub_residue, activity,
                 stmt, citation, evidence, annotations):
        super(ActivatingSubstitution, self).__init__(stmt, citation,
                                                     evidence, annotations)
        self.monomer_name = monomer_name
        self.wt_residue = wt_residue
        self.pos = pos
        self.sub_residue = sub_residue
        self.activity = activity

    def __str__(self):
        return ("ActivatingSubstitution(%s, %s, %s, %s, %s)" %
                (self.monomer_name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity))

class RasGef(Statement):
    """Statement representing the activation of a GTP-bound protein
    upon Gef activity."""

    def __init__(self, gef_name, gef_activity, ras_name,
                 stmt, citation, evidence, annotations):
        super(RasGef, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gef_name = gef_name
        self.gef_activity = gef_activity
        self.ras_name = ras_name

    def monomers(self, model, agent_set):
        gef = agent_set.get_create_agent(self.gef_name)
        gef.create_site(self.gef_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras_name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble(self, model, agent_set):
        try:
            kf_gef = model.parameters['kf_gef']
        except KeyError:
            kf_gef = Parameter('kf_gef', 1e-6)
            model.add_component(kf_gef)

        gef = model.monomers[self.gef_name]
        ras = model.monomers[self.ras_name]

        r = Rule('%s_activates_%s' %
                 (self.gef_name, self.ras_name),
                 gef(**{self.gef_activity: 'active'}) +
                 ras(**{'GtpBound': 'inactive'}) >>
                 gef(**{self.gef_activity: 'active'}) +
                 ras(**{'GtpBound': 'active'}),
                 kf_gef)
        model.add_component(r)

    def __str__(self):
        return ("RasGef(%s, %s, %s)" %
                (self.gef_name, self.gef_activity, self.ras_name))


class RasGap(Statement):
    """Statement representing the inactivation of a GTP-bound protein
    upon Gap activity."""
    def __init__(self, gap_name, gap_activity, ras_name,
                 stmt, citation, evidence, annotations):
        super(RasGap, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gap_name = gap_name
        self.gap_activity = gap_activity
        self.ras_name = ras_name

    def monomers(self, model, agent_set):
        gap = agent_set.get_create_agent(self.gap_name)
        gap.create_site(self.gap_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras_name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble(self, model, agent_set):
        try:
            kf_gap = model.parameters['kf_gap']
        except KeyError:
            kf_gap = Parameter('kf_gap', 1e-6)
            model.add_component(kf_gap)

        gap = model.monomers[self.gap_name]
        ras = model.monomers[self.ras_name]

        r = Rule('%s_inactivates_%s' %
                 (self.gap_name, self.ras_name),
                 gap(**{self.gap_activity: 'active'}) +
                 ras(**{'GtpBound': 'active'}) >>
                 gap(**{self.gap_activity: 'active'}) +
                 ras(**{'GtpBound': 'inactive'}),
                 kf_gap)
        model.add_component(r)

    def __str__(self):
        return ("RasGap(%s, %s, %s)" %
                (self.gap_name, self.gap_activity, self.ras_name))


class Complex(Statement):
    """Statement representing complex formation between a set of members"""
    def __init__(self, members):
        self.members = members

    def monomers(self, model, agent_set):
        """In this (very simple) implementation, proteins in a complex are
        each given site names corresponding to each of the other members
        of the complex. So the resulting complex is "fully connected" in
        that each is specified as bound to all the others."""
        for gene_name in self.members:
            gene_mono = agent_set.get_create_agent(gene_name)
            # Specify a binding site for each of the other complex members
            # bp = abbreviation for "binding partner"
            for bp_name in self.members:
                # The protein doesn't bind to itself!
                if gene_name == bp_name:
                    continue
                gene_mono.create_site(bp_name)

    def assemble(self, model, agent_set):
        # Get the rate parameter
        try:
            kf_bind = model.parameters['kf_bind']
        except KeyError:
            kf_bind = Parameter('kf_bind', 1e-6)
            model.add_component(kf_bind)

        # Make a rule name
        rule_name = '_'.join(self.members)
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
        for gene_name in self.members:
            mono = model.monomers[gene_name]
            # Specify free and bound states for binding sites for each of
            # the other complex members
            # (bp = abbreviation for "binding partner")
            free_site_dict = {}
            bound_site_dict = {}
            for bp_name in self.members:
                # The protein doesn't bind to itself!
                if gene_name == bp_name:
                    continue
                # Check to see if we've already created a bond index for these
                # two binding partners
                bp_set = ImmutableSet([gene_name, bp_name])
                if bp_set in bond_indices:
                    bond_ix = bond_indices[bp_set]
                # If we haven't see this pair of proteins yet, add a new bond
                # index to the dict
                else:
                    bond_ix = bond_counter
                    bond_indices[bp_set] = bond_ix
                    bond_counter += 1
                # Fill in the entries for the site dicts
                free_site_dict[bp_name] = None
                bound_site_dict[bp_name] = bond_ix
            # Build up the left- and right-hand sides of the rule from
            # monomer patterns with the appropriate site dicts
            mp_free = mono(**free_site_dict)
            mp_bound = mono(**bound_site_dict)
            lhs = lhs + mp_free
            rhs = rhs % mp_bound
        # Finally, create the rule and add it to the model
        try:
            rule = Rule(rule_name, lhs <> rhs, kf_bind, kf_bind)
            model.add_component(rule)
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

    def __str__(self):
        return ("Complex(%s)" % self.members)
