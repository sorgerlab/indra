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

def get_create_parameter(model, name, value, unique=True):
    """Return parameter with given name, creating it if needed.

    If unique is false and the parameter exists, the value is not changed; 
    if it does not exist, it will be created. If unique is true then upon conflic 
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

class UnknownPolicyException(Exception):
    pass


class Agent(object):
    def __init__(self, name, mods=None, mod_sites=None, 
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
        if db_refs is None:
            self.db_refs = {}
        else:
            self.db_refs = db_refs
        self.bound_to = bound_to
        self.bound_neg = bound_neg

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

class Statement(object):
    """The parent class of all statements"""
    def __init__(self, stmt, citation, evidence, annotations):
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
    def __init__(self, enz, mod, mod_pos, stmt,
                 citation, evidence, annotations):
        super(SelfModification, self).__init__(stmt, citation, evidence,
                                           annotations)
        self.enz = enz
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        return ("%s(%s, %s, %s)" %
                (type(self).__name__, self.enz.name, self.mod, self.mod_pos))

    def __eq__(self, other):
        if self.enz == other.enz and\
            self.mod == other.mod and\
            self.mod_pos == other.mod_pos:
            return True
        else:
            return False

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
            enz.create_site(self.enz.bound_to)

    def assemble_interactions_only(self, model, agent_set):
        kf_bind = get_create_parameter(model, 'kf_bind', 1.0, unique=False)
        kr_bind = get_create_parameter(model, 'kr_bind', 1.0, unique=False)

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
                     kf_bind, kr_bind)
            model.add_component(r)
        # If this rule is already in the model, issue a warning and continue
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

    def assemble_one_step(self, model, agent_set):
        enz = model.monomers[self.enz.name]
        sub = model.monomers[self.sub.name]

        param_name = 'kf_' + enz.name[0].lower() + sub.name[0].lower() + '_phos'
        kf_phospho = get_create_parameter(model, param_name, 1.0)

        # See NOTE in monomers_one_step
        site = site_name(self)[0]

        rule_name = '%s_phospho_%s_%s' % (self.enz.name, self.sub.name, site)
        try:
            # Iterate over all of the activating modification states of the
            # kinase
            enz_act_mods = agent_set[self.enz.name].activating_mods
            if enz_act_mods:
                for act_mod_pattern in [a.copy() for a in enz_act_mods]:
                    # Here we make the assumption that the binding site
                    # is simply named after the binding partner
                    if self.enz.bound_to:
                        if self.enz.bound_neg:
                            act_mod_pattern[self.enz.bound_to] = None
                            enz_pattern = enz(**act_mod_pattern)
                        else:
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
                    if self.enz.bound_neg:
                        enz_pattern = enz(**{self.enz.bound_to:None})
                    else:
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

# Autophosphorylation happens when a protein phosphorylates itself.
# A more precise name for this is cis-autophosphorylation.
class Autophosphorylation(SelfModification):
    def monomers_interactions_only(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        enz.create_site(site_name(self)[0], ('u', 'p'))

    def monomers_one_step(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        # NOTE: This assumes that a Phosphorylation statement will only ever
        # involve a single phosphorylation site on the substrate (typically
        # if there is more than one site, they will be parsed into separate
        # Phosphorylation statements, i.e., phosphorylation is assumed to be
        # distributive. If this is not the case, this assumption will need to
        # be revisited.
        enz.create_site(site_name(self)[0], ('u', 'p'))

        if self.enz.bound_to:
            enz_bound = agent_set.get_create_agent(self.enz.bound_to)
            enz_bound.create_site(self.enz.name)
            enz.create_site(self.enz.bound_to)

    def assemble_interactions_only(self, model, agent_set):
        self.assemble_one_step(model, agent_set)

    def assemble_one_step(self, model, agent_set):
        enz = model.monomers[self.enz.name]
        param_name = 'kf_' + enz.name[0].lower() + '_autophos'
        kf_autophospho = get_create_parameter(model, param_name, 1e-3)


        # See NOTE in monomers_one_step
        site = site_name(self)[0]

        rule_name = '%s_autophospho_%s_%s' % (self.enz.name, self.enz.name, site)
        try:
            # Iterate over all of the activating modification states of the
            # kinase
            enz_act_mods = agent_set[self.enz.name].activating_mods
            if enz_act_mods:
                for mod_pattern in [a.copy() for a in enz_act_mods]:
                    # Here we make the assumption that the binding site
                    # is simply named after the binding partner
                    left_pattern = mod_pattern.copy()
                    right_pattern = mod_pattern.copy()
                    left_pattern[site] = 'u'
                    right_pattern[site] = 'p'
                    if self.enz.bound_to:
                        if self.enz.bound_neg:
                            left_pattern[self.enz.bound_to] = None
                            right_pattern[self.enz.bound_to] = None
                            left_enz = enz(**left_pattern)
                            right_enz = enz(**right_pattern)
                        else:
                            left_pattern[self.enz.bound_to] = 1
                            right_pattern[self.enz.bound_to] = 1
                            enz_bound = model.monomers[self.enz.bound_to]
                            left_enz = enz(**left_pattern) % \
                                           enz_bound(**{self.enz.name:1})
                            right_enz = enz(**right_pattern) % \
                                            enz_bound(**{self.enz.name:1})
                    else:
                        left_enz = enz(**left_pattern)
                        right_enz = enz(**right_pattern)

                    r = Rule(rule_name, left_enz >> right_enz, kf_autophospho)
                    model.add_component(r)
            # If there are no known activity modifications, we take this
            # statement as given and allow the enzyme to phosphorylate itself
            # unconditionally
            else:
                left_pattern = {site: 'u'}
                right_pattern = {site: 'p'}
                if self.enz.bound_to:
                    if self.enz.bound_neg:
                        enz_bound = model.monomers[self.enz.bound_to]
                        left_pattern[self.enz.bound_to] = None
                        right_pattern[self.enz.bound_to] = None
                        left_enz = enz(**left_pattern)
                        right_enz = enz(**right_pattern)
                    else:
                        enz_bound = model.monomers[self.enz.bound_to]
                        left_pattern[self.enz.bound_to] = 1
                        right_pattern[self.enz.bound_to] = 1
                        left_enz = enz(**left_pattern) % \
                                    enz_bound(**{self.enz.name:1})
                        right_enz = enz(**right_pattern) % \
                                    enz_bound(**{self.enz.name:1})
                else:
                    left_enz = enz(**left_pattern)
                    right_enz = enz(**right_pattern)
                r = Rule(rule_name, left_enz >> right_enz, kf_autophospho)
                model.add_component(r)

        # If this rule is already in the model, issue a warning and continue
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

# Transphosphorylation assumes that a kinase is already bound to 
# a substrate (usually of the same molecular species), and phosphorylates 
# it in an intra-molecular fashion. The enz property of the statement must 
# have exactly one bound_to property, and we assume that enz phosphorylates 
# this bound_to molecule. The bound_neg property is ignored here.
class Transphosphorylation(SelfModification):
    def monomers_interactions_only(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        # Assume there is exactly one bound_to species
        sub = agent_set.get_create_agent(self.enz.bound_to)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def monomers_one_step(self, agent_set):
        enz = agent_set.get_create_agent(self.enz.name)
        # NOTE: This assumes that a Phosphorylation statement will only ever
        # involve a single phosphorylation site on the substrate (typically
        # if there is more than one site, they will be parsed into separate
        # Phosphorylation statements, i.e., phosphorylation is assumed to be
        # distributive. If this is not the case, this assumption will need to
        # be revisited. 
        # Assume there is exactly one bound_to species
        sub = agent_set.get_create_agent(self.enz.bound_to)
        sub.create_site(site_name(self)[0], ('u', 'p'))

        if self.enz.bound_to:
            enz_bound = agent_set.get_create_agent(self.enz.bound_to)
            enz_bound.create_site(self.enz.name)
            enz.create_site(self.enz.bound_to)

    def assemble_interactions_only(self, model, agent_set):
        self.assemble_one_step(model, agent_set)

    def assemble_one_step(self, model, agent_set):
        enz = model.monomers[self.enz.name]
        sub = model.monomers[self.enz.bound_to]
        param_name = 'kf_' + enz.name[0].lower() + sub.name[0].lower() + '_transphos'
        kf  = get_create_parameter(model, param_name, 1e-3)
        
        # See NOTE in monomers_one_step
        site = site_name(self)[0]

        rule_name = '%s_transphospho_%s_%s' % (self.enz.name, self.enz.bound_to, site)
        try:
            # Iterate over all of the activating modification states of the
            # kinase
            enz_act_mods = agent_set[self.enz.name].activating_mods
            if enz_act_mods:
                for mod_pattern in [a.copy() for a in enz_act_mods]:
                    # Here we make the assumption that the binding site
                    # is simply named after the binding partner
                    left_pattern = mod_pattern.copy()
                    right_pattern = mod_pattern.copy()
                    left_pattern[self.enz.bound_to] = 1
                    right_pattern[self.enz.bound_to] = 1
                    enz_bound = model.monomers[self.enz.bound_to]
                    left_enz = enz(**left_pattern) % \
                                enz_bound(**{self.enz.name: 1, site: 'u'})
                    right_enz = enz(**right_pattern) % \
                                enz_bound(**{self.enz.name: 1, site: 'p'})

                    r = Rule(rule_name, left_enz >> right_enz, kf_autophospho)
                    model.add_component(r)
            # If there are no known activity modifications, we take this
            # statement as given and allow the enzyme to phosphorylate itself
            # unconditionally
            else:
                left_pattern = {}
                right_pattern = {}
                enz_bound = model.monomers[self.enz.bound_to]
                left_pattern[self.enz.bound_to] = 1
                right_pattern[self.enz.bound_to] = 1
                left_enz = enz(**left_pattern) % \
                            enz_bound(**{self.enz.name: 1, site: 'u'})
                right_enz = enz(**right_pattern) % \
                            enz_bound(**{self.enz.name: 1, site: 'p'})
                r = Rule(rule_name, left_enz >> right_enz, kf_autophospho)
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
        subj = agent_set.get_create_agent(self.subj.name)
        subj.create_site(active_site_names[self.subj_activity])
        obj = agent_set.get_create_agent(self.obj.name)
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
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        subj = agent_set.get_create_agent(self.subj.name)
        subj.create_site(self.subj_activity, ('inactive', 'active'))
        obj = agent_set.get_create_agent(self.obj.name)
        obj.create_site(self.obj_activity, ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        subj = model.monomers[self.subj.name]
        obj = model.monomers[self.obj.name]

        param_name = 'kf_' + subj.name[0].lower() + obj.name[0].lower() + '_act'
        kf_one_step_activate = \
                       get_create_parameter(model, param_name, 1e-6)

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
        phos = agent_set.get_create_agent(self.phos.name)
        phos.create_site(active_site_names['Phosphatase'])
        sub = agent_set.get_create_agent(self.sub.name)
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
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        phos = agent_set.get_create_agent(self.phos.name)
        sub = agent_set.get_create_agent(self.sub.name)
        sub.create_site(site_name(self)[0], ('u', 'p'))

    def assemble_one_step(self, model, agent_set):
        phos = model.monomers[self.phos.name]
        sub = model.monomers[self.sub.name]

        param_name = 'kf_' + phos.name[0].lower() + sub.name[0].lower() + '_dephos'
        kf_dephospho = get_create_parameter(model, param_name, 1e-6)
        
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
        agent = agent_set.get_create_agent(self.monomer.name)
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
                 stmt, citation, evidence, annotations):
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
                 stmt, citation, evidence, annotations):
        super(RasGef, self).__init__(stmt, citation, evidence,
                                     annotations)
        self.gef = gef
        self.gef_activity = gef_activity
        self.ras = ras

    def __eq__(self, other):
        if self.gef == other.gef and\
            self.gef_activity == other.gef_activity and\
            self.ras == other.ras:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        gef = agent_set.get_create_agent(self.gef.name)
        gef.create_site('gef_site')
        ras = agent_set.get_create_agent(self.ras.name)
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
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        gef = agent_set.get_create_agent(self.gef.name)
        gef.create_site(self.gef_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        gef = model.monomers[self.gef.name]
        ras = model.monomers[self.ras.name]

        param_name = 'kf_' + gef.name[0].lower() + ras.name[0].lower() + '_gef'
        kf_gef = get_create_parameter(model, param_name, 1e-6)

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
    
    def __eq__(self, other):
        if isinstance(other, RasGap) and\
            self.gap == other.gap and\
            self.gap_activity == other.gap_activity and\
            self.ras == other.ras:
            return True
        else:
            return False

    def monomers_interactions_only(self, agent_set):
        gap = agent_set.get_create_agent(self.gap.name)
        gap.create_site('gap_site')
        ras = agent_set.get_create_agent(self.ras.name)
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
        model.add_component(r)

    def monomers_one_step(self, agent_set):
        gap = agent_set.get_create_agent(self.gap.name)
        gap.create_site(self.gap_activity, ('inactive', 'active'))
        ras = agent_set.get_create_agent(self.ras.name)
        ras.create_site('GtpBound', ('inactive', 'active'))

    def assemble_one_step(self, model, agent_set):
        gap = model.monomers[self.gap.name]
        ras = model.monomers[self.ras.name]
        
        param_name = 'kf_' + gap.name[0].lower() + ras.name[0].lower() + '_gap'
        kf_gap = get_create_parameter(model, param_name, 1e-6)

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

    def __eq__(self, other):
        # TODO: find equality for different orders of members too
        if not isinstance(other, Complex):
            return False
        for (m1, m2) in zip(self.members, other.members):
            if not m1 == m2:
                return False
        return True

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
            # Add sites for agent modifications
            # TODO: This assumes phosphorylation, but in principle
            # it could be some other modification
            for m, mp in zip(member.mods, member.mod_sites):
                mod = abbrevs[m]
                mod_pos = mp if mp is not None else ''
                mod_site = ('%s%s' % (mod, mod_pos))
                gene_mono.create_site(mod_site, ['u', 'p'])
            # Add binding sites when agent is bound to something
            # This assumes that there is only one binding partner
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
        abbr_name = ''.join([m.name[0].lower() for m in self.members])
        kf_bind = get_create_parameter(model, 'kf_' + abbr_name + '_bind', 1e-6)
        kr_bind = get_create_parameter(model, 'kr_' + abbr_name + '_bind', 1e-6)

        # Make a rule name
        name_components = []
        for m in self.members:
            if m.bound_to:
                if m.bound_neg:
                    name_components.append(m.name + '_n' + m.bound_to)
                else:
                    name_components.append(m.name + '_' + m.bound_to)
            else:
                name_components.append(m.name)
        rule_name = '_'.join(name_components) + '_bind'
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
            
            # Add the pattern for the modifications of the member
            # TODO: This is specific to phosphorylation but we should be 
            # able to support other types as well
            for m, mp in zip(member.mods, member.mod_sites):
                mod = abbrevs[m]
                mod_pos = mp if mp is not None else ''
                mod_site = ('%s%s' % (mod, mod_pos)) 
                left_site_dict[mod_site] = 'p'
                right_site_dict[mod_site] = 'p'

            # Add the pattern for the member being bound
            if member.bound_to:
                if member.bound_neg:
                    bound_name = member.bound_to
                    left_site_dict[bound_name] = None
                    right_site_dict[bound_name] = None
                    left_pattern = mono(**left_site_dict)
                    right_pattern = mono(**right_site_dict)
                else:
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
            rule = Rule(rule_name, lhs <> rhs, kf_bind, kr_bind)
            model.add_component(rule)
        except ComponentDuplicateNameError:
            msg = "Rule %s already in model! Skipping." % rule_name
            warnings.warn(msg)

    def __str__(self):
        return ("Complex(%s)" % [m.name for m in self.members])
