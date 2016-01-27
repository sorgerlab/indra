import warnings
from pysb import *
from pysb import ReactionPattern, ComplexPattern, ComponentDuplicateNameError


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

def add_rule_to_model(model, rule):
    try:
        model.add_component(rule)
    # If this rule is already in the model, issue a warning and continue
    except ComponentDuplicateNameError:
        msg = "Rule %s already in model! Skipping." % rule.name
        warnings.warn(msg)

def get_activating_mods(agent, agent_set):
    act_mods = agent_set[agent.name].activating_mods
    if not act_mods:
        act_mods = [{}]
    return act_mods

def get_binding_site_name(name):
    binding_site = name.lower()
    return binding_site


class Agent(object):
    def __init__(self, name, mods=None, mod_sites=None, active=None,
                 bound_to=None, bound_neg=None, db_refs=None):
        self.name = name
        if mods is None:
            self.mods = []
        else:
            self.mods = mods
        if mod_sites is None:
            self.mod_sites = []
        else:
            self.mod_sites = mod_sites
        self.bound_to = bound_to
        self.bound_neg = bound_neg
        self.active = active
        if db_refs is None:
            self.db_refs = {}
        else:
            self.db_refs = db_refs

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        attr_strs = []
        if self.mods:
            attr_strs.append('mods: %s' % self.mods)
        if self.mod_sites:
            attr_strs.append('mod_sites: %s' % self.mod_sites)
        if self.active:
            attr_strs.append('active: %s' % self.active)
        if self.bound_to:
            attr_strs.append('bound_to: %s' % self.bound_to)
        if self.bound_neg:
            attr_strs.append('bound_neg: %s' % self.bound_neg)
        if self.db_refs:
            attr_strs.append('db_refs: %s' % self.db_refs)
        attr_str = ', '.join(attr_strs)
        return '%s(%s)' % (self.name, attr_str)

class Statement(object):
    """The parent class of all statements"""
    def __init__(self, stmt=None, citation=None, evidence=None, 
                 annotations=None):
        self.stmt = stmt
        self.citation = citation
        self.evidence = evidence
        self.annotations = annotations

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self.citation == other.citation and\
            self.evidence == other.evidence and\
            self.annotations == other.annotations and\
            self.stmt == other.stmt:
            return True
        else:
            return False


class Modification(Statement):
    """Generic statement representing the modification of a protein"""
    def __init__(self, enz, sub, mod, mod_pos, stmt=None,
                 citation=None, evidence=None, annotations=None):
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

    def __eq__(self, other):
        if isinstance(other, Modification) and\
            self.enz == other.enz and\
            self.sub == other.sub and\
            self.mod == other.mod and\
            self.mod_pos == other.mod_pos:
            return True
        else:
            return False

class SelfModification(Statement):
    """Generic statement representing the self modification of a protein"""
    def __init__(self, enz, mod, mod_pos, stmt=None,
                 citation=None, evidence=None, annotations=None):
        super(SelfModification, self).__init__(stmt, citation, evidence,
                                           annotations)
        self.enz = enz
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        return ("%s(%s, %s, %s)" %
                (type(self).__name__, self.enz.name, self.mod, self.mod_pos))

    def __eq__(self, other):
        if isinstance(other, SelfModification) and\
            self.enz == other.enz and\
            self.mod == other.mod and\
            self.mod_pos == other.mod_pos:
            return True
        else:
            return False


class Phosphorylation(Modification):
    """Phosphorylation modification"""

# Autophosphorylation happens when a protein phosphorylates itself.
# A more precise name for this is cis-autophosphorylation.
class Autophosphorylation(SelfModification):
    def monomers_interactions_only(self, agent_set):
        enz = agent_set.get_create_base_agent(self.enz)
        enz.create_site(site_name(self)[0], ('u', 'p'))

    def monomers_one_step(self, agent_set):
        enz = agent_set.get_create_base_agent(self.enz)
        # NOTE: This assumes that a Phosphorylation statement will only ever
        # involve a single phosphorylation site on the substrate (typically
        # if there is more than one site, they will be parsed into separate
        # Phosphorylation statements, i.e., phosphorylation is assumed to be
        # distributive. If this is not the case, this assumption will need to
        # be revisited.
        enz.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_interactions_only(self, model, agent_set):
        self.assemble_one_step(model, agent_set)

    def assemble_one_step(self, model, agent_set):
        param_name = 'kf_' + self.enz.name[0].lower() + '_autophos'
        kf_autophospho = get_create_parameter(model, param_name, 1e-3)


        # See NOTE in monomers_one_step
        site = site_name(self)[0]
        pattern_unphos = get_complex_pattern(model, self.enz, agent_set, extra_fields={site: 'u'})
        pattern_phos = get_complex_pattern(model, self.enz, agent_set, extra_fields={site: 'p'})

        rule_name = '%s_autophospho_%s_%s' % (self.enz.name, self.enz.name, site)
        r = Rule(rule_name, pattern_unphos >> pattern_phos, kf_autophospho)
        add_rule_to_model(model, r)

# Transphosphorylation assumes that a kinase is already bound to 
# a substrate (usually of the same molecular species), and phosphorylates 
# it in an intra-molecular fashion. The enz property of the statement must 
# have exactly one bound_to property, and we assume that enz phosphorylates 
# this bound_to molecule. The bound_neg property is ignored here.
class Transphosphorylation(SelfModification):
    def monomers_interactions_only(self, agent_set):
        enz = agent_set.get_create_base_agent(self.enz)
        # Assume there is exactly one bound_to species
        sub = agent_set.get_create_base_agent(self.enz)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def monomers_one_step(self, agent_set):
        enz = agent_set.get_create_base_agent(self.enz)
        # NOTE: This assumes that a Phosphorylation statement will only ever
        # involve a single phosphorylation site on the substrate (typically
        # if there is more than one site, they will be parsed into separate
        # Phosphorylation statements, i.e., phosphorylation is assumed to be
        # distributive. If this is not the case, this assumption will need to
        # be revisited. 
        sub = agent_set.get_create_base_agent(Agent(self.enz.bound_to))
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_interactions_only(self, model, agent_set):
        self.assemble_one_step(model, agent_set)

    def assemble_one_step(self, model, agent_set):
        param_name = 'kf_' + self.enz.name[0].lower() + self.enz.bound_to[0].lower() + '_transphos'
        kf  = get_create_parameter(model, param_name, 1e-3)
        
        site = site_name(self)[0]
        enz_pattern = get_complex_pattern(model, self.enz, agent_set)
        sub_unphos = get_complex_pattern(model, Agent(self.enz.bound_to), 
            agent_set, extra_fields = {site: 'u'})
        sub_phos = get_complex_pattern(model, Agent(self.enz.bound_to), 
            agent_set, extra_fields = {site: 'p'})
        
        rule_name = '%s_transphospho_%s_%s' % (self.enz.name, self.enz.bound_to, site)    
        r = Rule(rule_name, enz_pattern % sub_unphos >>\
                        enz_pattern % sub_phos, kf)
        add_rule_to_model(model, r)

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
                 obj_activity, stmt=None, citation=None, evidence=None,
                 annotations=None):
        super(ActivityActivity, self).__init__(stmt,
                                               citation, evidence, annotations)
        self.subj = subj
        self.subj_activity = subj_activity
        self.obj = obj
        self.obj_activity = obj_activity
        self.relationship = relationship

    def __eq__(self, other):
        if isinstance(other, ActivityActivity) and\
            self.subj == other.subj and\
            self.subj_activity == other.subj_activity and\
            self.obj == other.obj and\
            self.obj_activity == other.obj_activity:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        subj = agent_set.get_create_base_agent(self.subj)
        subj.create_site(active_site_names[self.subj_activity])
        obj = agent_set.get_create_base_agent(self.obj)
        obj.create_site(active_site_names[self.obj_activity])
        obj.create_site(default_mod_site_names[self.subj_activity])

    def assemble_interactions_only(self, model):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
        subj = model.monomers[self.subj.name]
        obj = model.monomers[self.obj.name]
        subj_active_site = active_site_names[self.subj_activity]
        obj_mod_site = default_mod_site_names[self.subj_activity]
        r = Rule('%s_%s_activates_%s_%s' %
                 (self.subj.name, self.subj_activity, self.obj.name,
                  self.obj_activity),
                 subj(**{subj_active_site: None}) +
                 obj(**{obj_mod_site: None}) >>
                 subj(**{subj_active_site: 1}) %
                 obj(**{obj_mod_site: 1}),
                 kf_bind)
        add_rule_to_model(model, r)

    def monomers_one_step(self, agent_set):
        subj = agent_set.get_create_base_agent(self.subj)
        subj.create_site(self.subj_activity, ('inactive', 'active'))
        obj = agent_set.get_create_base_agent(self.obj)
        obj.create_site(self.obj_activity, ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        subj_pattern = get_complex_pattern(model, self.subj, agent_set, 
            extra_fields={self.subj_activity: 'active'})
        obj_inactive = get_complex_pattern(model, self.obj, agent_set, 
            extra_fields={self.obj_activity: 'inactive'})
        obj_active = get_complex_pattern(model, self.obj, agent_set, 
            extra_fields={self.obj_activity: 'active'})
        
        param_name = 'kf_' + self.subj.name[0].lower() +\
                            self.obj.name[0].lower() + '_act'
        kf_one_step_activate = \
                       get_create_parameter(model, param_name, 1e-6)

        rule_name = '%s_%s_activates_%s_%s' %\
            (self.subj.name, self.subj_activity, self.obj.name, self.obj_activity)
        if self.relationship == 'increases':
           r = Rule(rule_name,                
                subj_pattern + obj_inactive >> subj_pattern + obj_active,
                kf_one_step_activate)
        else:
           r = Rule(rule_name,                
                subj_pattern + obj_active >> subj_pattern + obj_inactive,
                kf_one_step_activate)

        add_rule_to_model(model, r)

    def __str__(self):
        return ("%s(%s, %s, %s, %s, %s)" %
                (type(self).__name__, self.subj.name, self.subj_activity,
                 self.relationship, self.obj.name, self.obj_activity))

class RasGtpActivityActivity(ActivityActivity):
    pass

class Dephosphorylation(Statement):
    def __init__(self, phos, sub, mod, mod_pos, stmt=None,
                 citation=None, evidence=None, annotations=None):
        super(Dephosphorylation, self).__init__(stmt, citation,
                                                evidence, annotations)
        self.phos = phos
        self.sub = sub
        self.mod = mod
        self.mod_pos = mod_pos
    
    def __eq__(self, other):
        if isinstance(other, Dephosphorylation) and\
            self.phos == other.phos and\
            self.sub == other.sub and\
            self.mod == other.mod and\
            self.mod_pos == other.mod_pos:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        phos = agent_set.get_create_base_agent(self.phos)
        phos.create_site(active_site_names['Phosphatase'])
        sub = agent_set.get_create_base_agent(self.sub)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
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
        add_rule_to_model(model, r)

    def monomers_one_step(self, agent_set):
        phos = agent_set.get_create_base_agent(self.phos)
        sub = agent_set.get_create_base_agent(self.sub)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_one_step(self, model, agent_set):
        param_name = 'kf_' + self.phos.name[0].lower() +\
                    self.sub.name[0].lower() + '_dephos'
        kf_dephospho = get_create_parameter(model, param_name, 1e-6)
        
        site = site_name(self)[0]
        phos_pattern = get_complex_pattern(model, self.phos, agent_set)
        sub_phos = get_complex_pattern(model, self.sub, agent_set, 
            extra_fields={site: 'p'})
        sub_unphos = get_complex_pattern(model, self.sub, agent_set, 
            extra_fields={site: 'u'})

        r = Rule('%s_dephospho_%s_%s' %
                 (self.phos.name, self.sub.name, site),
                 phos_pattern + sub_phos >>
                 phos_pattern + sub_unphos,
                 kf_dephospho)
        add_rule_to_model(model, r)

    def monomers_two_step(self, agent_set):
        phos = agent_set.get_create_base_agent(self.phos)
        sub = agent_set.get_create_base_agent(self.sub)
        sub.create_site(site_name(self)[0], ('u', 'p'))

        # Create site for binding the substrate
        phos.create_site(get_binding_site_name(sub.name))
        sub.create_site(get_binding_site_name(phos.name))

    def assemble_two_step(self, model, agent_set):
        sub_bs = get_binding_site_name(self.sub.name)
        phos_bs = get_binding_site_name(self.phos.name)
        phos_bound = get_complex_pattern(model, self.phos, agent_set,
            extra_fields = {sub_bs: 1})
        phos_unbound = get_complex_pattern(model, self.phos, agent_set,
            extra_fields = {sub_bs: None})
        sub_pattern = get_complex_pattern(model, self.sub, agent_set)
        
        param_name = 'kf_' + self.phos.name[0].lower() +\
            self.sub.name[0].lower() + '_bind'
        kf_bind = get_create_parameter(model, param_name, 1e-6)
        param_name = 'kr_' + self.phos.name[0].lower() +\
            self.sub.name[0].lower() + '_bind'
        kr_bind = get_create_parameter(model, param_name, 1e-3)
        param_name = 'kc_' + self.phos.name[0].lower() +\
            self.sub.name[0].lower() + '_dephos'
        kf_phospho = get_create_parameter(model, param_name, 1e-3)
        
        site = site_name(self)[0]

        phos_act_mods = get_activating_mods(self.phos, agent_set)
        for i, am in enumerate(phos_act_mods):
            rule_name = '%s_dephos_bind_%s_%s_%d' %\
                (self.phos.name, self.sub.name, site, i+1)
            r = Rule(rule_name,
                phos_unbound(am) +\
                sub_pattern(**{site: 'p', phos_bs: None}) >>
                phos_bound(am) %\
                sub_pattern(**{site: 'p', phos_bs: 1}),
                kf_bind, kr_bind)
            add_rule_to_model(model, r)
        
            rule_name = '%s_dephos_%s_%s_%d' %\
                (self.phos.name, self.sub.name, site, i+1)
            r = Rule(rule_name,
                phos_bound(am) %\
                    sub_pattern(**{site: 'p', phos_bs: 1}) >>
                phos_unbound(am) +\
                    sub_pattern(**{site: 'u', phos_bs: None}),
                kf_phospho)
            add_rule_to_model(model, r)
        
        rule_name = '%s_dissoc_%s' % (self.phos.name, self.sub.name)
        r = Rule(rule_name, model.monomers[self.phos.name](**{sub_bs: 1}) %\
                 model.monomers[self.sub.name](**{phos_bs: 1}) >>
                 model.monomers[self.phos.name](**{sub_bs: None}) +\
                 model.monomers[self.sub.name](**{phos_bs: None}), kr_bind)
        add_rule_to_model(model, r)

    def __str__(self):
        return ("Dephosphorylation(%s, %s, %s, %s)" %
                (self.phos.name, self.sub.name, self.mod, self.mod_pos))

class ActivityModification(Statement):
    """Statement representing the activation of a protein as a result
    of a residue modification"""
    def __init__(self, monomer, mod, mod_pos, relationship, activity,
                 stmt=None, citation=None, evidence=None, annotations=None):
        super(ActivityModification, self).__init__(stmt, citation,
                                                   evidence, annotations)
        self.monomer = monomer
        self.mod = mod
        self.mod_pos = mod_pos
        self.relationship = relationship
        self.activity = activity
    
    def __eq__(self, other):
        if isinstance(other, ActivityModification) and\
            self.monomer == other.monomer and\
            self.mod == other.mod and\
            self.mod_pos == other.mod_pos and\
            self.relationship == other.relationship and\
            self.activity == other.activity:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def monomers_one_step(self, agent_set):
        agent = agent_set.get_create_base_agent(self.monomer)
        sites = site_name(self)
        active_states = [states[m][1] for m in self.mod]

        activity_pattern = {}
        for i, s in enumerate(sites):
            site_states = states[self.mod[i]]
            active_state = site_states[1]
            agent.create_site(s, site_states)
            activity_pattern[s] = active_state

        # Add this activity modification explicitly to the agent's list
        # of activating modifications
        agent.add_activating_modification(activity_pattern)
        # Inactivating modifications will require a different treatment
        # of the resolution of when the agent is active
        if self.relationship == 'decreases':
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
    def __init__(self, monomer, wt_residue, pos, sub_residue, activity, rel,
                 stmt=None, citation=None, evidence=None, annotations=None):
        super(ActivatingSubstitution, self).__init__(stmt, citation,
                                                     evidence, annotations)
        self.monomer = monomer
        self.wt_residue = wt_residue
        self.pos = pos
        self.sub_residue = sub_residue
        self.activity = activity
        self.rel = rel
    
    def __eq__(self, other):
        if isinstance(other, ActivatingSubstitution) and\
            self.monomer == other.monomer and\
            self.wt_residue == other.wt_residue and\
            self.pos == other.pos and\
            self.sub_residue == other.sub_residue and\
            self.activity == other.activity:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        pass

    def assemble_interactions_only(self, model, agent_set):
        pass

    def __str__(self):
        return ("ActivatingSubstitution(%s, %s, %s, %s, %s, %s)" %
                (self.monomer.name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity, self.rel))

class RasGef(Statement):
    """Statement representing the activation of a GTP-bound protein
    upon Gef activity."""

    def __init__(self, gef, gef_activity, ras,
                 stmt=None, citation=None, evidence=None, annotations=None):
        super(RasGef, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gef = gef
        self.gef_activity = gef_activity
        self.ras = ras

    def __eq__(self, other):
        if isinstance(other, RasGef) and\
            self.gef == other.gef and\
            self.gef_activity == other.gef_activity and\
            self.ras == other.ras:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        gef = agent_set.get_create_base_agent(self.gef)
        gef.create_site('gef_site')
        ras = agent_set.get_create_base_agent(self.ras)
        ras.create_site('p_loop')

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
        gef = model.monomers[self.gef.name]
        ras = model.monomers[self.ras.name]
        r = Rule('%s_activates_%s' %
                 (self.gef.name, self.ras.name),
                 gef(**{'gef_site':None}) +
                 ras(**{'p_loop':None}) >>
                 gef(**{'gef_site': 1}) +
                 ras(**{'p_loop': 1}),
                 kf_bind)
        add_rule_to_model(model, r)

    def monomers_one_step(self, agent_set):
        gef = agent_set.get_create_base_agent(self.gef)
        gef.create_site(self.gef_activity, ('inactive', 'active'))
        ras = agent_set.get_create_base_agent(self.ras)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        gef_pattern = get_complex_pattern(model, self.gef, agent_set, 
            extra_fields={self.gef_activity: 'active'})
        ras_inactive = get_complex_pattern(model, self.ras, agent_set,
            extra_fields={'GtpBound': 'inactive'})
        ras_active = get_complex_pattern(model, self.ras, agent_set,
            extra_fields={'GtpBound': 'active'})

        param_name = 'kf_' + self.gef.name[0].lower() +\
                        self.ras.name[0].lower() + '_gef'
        kf_gef = get_create_parameter(model, param_name, 1e-6)

        r = Rule('%s_activates_%s' %
                 (self.gef.name, self.ras.name),
                 gef_pattern + ras_inactive >>
                 gef_pattern + ras_active,
                 kf_gef)
        add_rule_to_model(model, r)

    def __str__(self):
        return ("RasGef(%s, %s, %s)" %
                (self.gef.name, self.gef_activity, self.ras.name))


class RasGap(Statement):
    """Statement representing the inactivation of a GTP-bound protein
    upon Gap activity."""
    def __init__(self, gap, gap_activity, ras,
                 stmt=None, citation=None, evidence=None, annotations=None):
        super(RasGap, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gap = gap
        self.gap_activity = gap_activity
        self.ras = ras
    
    def __eq__(self, other):
        if isinstance(other, RasGap) and\
            self.gap == other.gap and\
            self.gap_activity == other.gap_activity and\
            self.ras == other.ras:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        gap = agent_set.get_create_base_agent(self.gap)
        gap.create_site('gap_site')
        ras = agent_set.get_create_base_agent(self.ras)
        ras.create_site('gtp_site')

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
        gap = model.monomers[self.gap.name]
        ras = model.monomers[self.ras.name]
        r = Rule('%s_inactivates_%s' %
                 (self.gap.name, self.ras.name),
                 gap(**{'gap_site': None}) +
                 ras(**{'gtp_site': None}) >>
                 gap(**{'gap_site': 1}) +
                 ras(**{'gtp_site': 1}),
                 kf_bind)
        add_rule_to_model(model, r)

    def monomers_one_step(self, agent_set):
        gap = agent_set.get_create_base_agent(self.gap)
        gap.create_site(self.gap_activity, ('inactive', 'active'))
        ras = agent_set.get_create_base_agent(self.ras)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        gap_pattern = get_complex_pattern(model, self.gap, agent_set, 
            extra_fields={self.gap_activity: 'active'})
        ras_inactive = get_complex_pattern(model, self.ras, agent_set,
            extra_fields={'GtpBound': 'inactive'})
        ras_active = get_complex_pattern(model, self.ras, agent_set,
            extra_fields={'GtpBound': 'active'})

        param_name = 'kf_' + self.gap.name[0].lower() +\
                        self.ras.name[0].lower() + '_gap'
        kf_gap = get_create_parameter(model, param_name, 1e-6)

        r = Rule('%s_deactivates_%s' %
                 (self.gap.name, self.ras.name),
                 gap_pattern + ras_active >>
                 gap_pattern + ras_inactive,
                 kf_gap)
        add_rule_to_model(model, r)

    def __str__(self):
        return ("RasGap(%s, %s, %s)" %
                (self.gap.name, self.gap_activity, self.ras.name))


class Complex(Statement):
    """Statement representing complex formation between a set of members"""
    def __init__(self, members, stmt=None, citation=None, 
                 evidence=None, annotations=None):
        super(Complex, self).__init__(stmt, citation, evidence, annotations)
        self.members = members

    def __eq__(self, other):
        # TODO: find equality for different orders of members too
        if not isinstance(other, Complex):
            return False
        for (m1, m2) in zip(self.members, other.members):
            if not m1 == m2:
                return False
        return True



    def __str__(self):
        return ("Complex(%s)" % [m.name for m in self.members])
