import warnings
from pysb import *
from rdf_to_pysb import abbrevs, states

def site_name(stmt):
    """Return a site name for a modification-type statement."""
    mod = abbrevs[stmt.mod]
    mod_pos = stmt.mod_pos if stmt.mod_pos is not None else ''
    return '%s%s' % (mod, mod_pos)

def get_create_monomer(model, name):
    """Return monomer with given name, creating it if needed."""
    monomer = model.monomers.get(name)
    if monomer is None:
        monomer = Monomer(name)
        model.add_component(monomer)
    return monomer

def create_site(monomer, site, states=None):
    if site not in monomer.sites:
        monomer.sites.append(site)
    if states is not None:
        monomer.site_states.setdefault(site, [])
        try:
            states = list(states)
        except TypeError:
            return
        add_site_states(monomer, site, states)

def add_site_states(monomer, site, states):
    for state in states:
        if state not in monomer.site_states[site]:
            monomer.site_states[site].append(state)

class Statement(object):
    def __init__(self, subj, obj, stmt, citation, evidence, annotations):
        self.subj = subj
        self.obj = obj
        self.stmt = stmt
        self.citation = citation
        self.evidence = evidence
        self.annotations = annotations

    def monomers(self, model):
        warnings.warn("%s.monomers not implemented" % self.__class__.__name__)

    def assemble(self, model):
        warnings.warn("%s.assemble not implemented" % self.__class__.__name__)

class Modification(Statement):
    def __init__(self, enz_name, sub_name, mod, mod_pos, subj, obj, stmt,
                 citation, evidence, annotations):
        super(Modification, self).__init__(subj, obj, stmt, citation, evidence,
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

    def monomers(self, model):
        enz = get_create_monomer(model, self.enz_name)
        create_site(enz, 'Kinase', ('inactive', 'active'))
        sub = get_create_monomer(model, self.sub_name)
        create_site(sub, site_name(self), ('u', 'p'))

    def assemble(self, model):
        try:
            kf_phospho = model.parameters['kf_phospho']
        except KeyError:
            kf_phospho = Parameter('kf_phospho', 1e-6)
            model.add_component(kf_phospho)

        enz = model.monomers[self.enz_name]
        sub = model.monomers[self.sub_name]

        site = site_name(self)
        r = Rule('%s_phospho_%s_%s' %
                 (self.enz_name, self.sub_name, site),
                 enz(Kinase='active') + sub(**{site:'u'}) >>
                 enz(Kinase='active') + sub(**{site:'p'}),
                 kf_phospho)
        model.add_component(r)

class Hydroxylation(Modification):
    pass

class Sumoylation(Modification):
    pass

class Acetylation(Modification):
    pass

class Ubiquitination(Modification):
    pass

class ActivityActivity(Statement):
    def __init__(self, subj_name, subj_activity, relationship, obj_name,
                 obj_activity, subj, obj, stmt, citation, evidence,
                 annotations):
        super(ActivityActivity, self).__init__(subj, obj, stmt,
                                               citation, evidence, annotations)
        self.subj_name = subj_name
        self.subj_activity = subj_activity
        self.obj_name = obj_name
        self.obj_activity = obj_activity
        self.relationship = relationship

    def monomers(self, model):
        subj = get_create_monomer(model, self.subj_name)
        create_site(subj, self.subj_activity, ('inactive','active'))
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
                 subj(**{self.subj_activity:'active'}) +
                 obj(**{self.obj_activity:'inactive'}) >>
                 subj(**{self.subj_activity:'active'}) +
                 obj(**{self.obj_activity:'active'}),
                 kf_one_step_activate)
        model.add_component(r)

    def __str__(self):
        return ("%s(%s, %s, %s, %s)" %
                (type(self).__name__, self.subj_name, self.subj_activity,
                 self.obj_name, self.obj_activity))

class Dephosphorylation(Statement):
    def __init__(self, phos_name, sub_name, mod, mod_pos, subj, obj, stmt,
                 citation, evidence, annotations):
        super(Dephosphorylation, self).__init__(subj, obj, stmt, citation,
                                                evidence, annotations)
        self.phos_name = phos_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def monomers(self, model):
        phos = get_create_monomer(model, self.phos_name)
        sub = get_create_monomer(model, self.sub_name)
        create_site(sub, site_name(self), ('u', 'p'))

    def assemble(self, model):
        try:
            kf_dephospho = model.parameters['kf_dephospho']
        except KeyError:
            kf_dephospho = Parameter('kf_dephospho', 1e-6)
            model.add_component(kf_dephospho)

        phos = model.monomers[self.phos_name]
        sub = model.monomers[self.sub_name]

        site = site_name(self)
        r = Rule('%s_dephospho_%s_%s' %
                 (self.phos_name, self.sub_name, site),
                 phos() + sub(**{site:'p'}) >>
                 phos() + sub(**{site:'u'}),
                 kf_dephospho)
        model.add_component(r)


    def __str__(self):
        return ("Dehosphorylation(%s, %s, %s, %s)" %
                (self.phos_name, self.sub_name, self.mod, self.mod_pos))

class ActivityModification(Statement):
    def __init__(self, monomer_name, mod, mod_pos, relationship, activity,
                 subj, obj, stmt, citation, evidence, annotations):
        super(ActivityModification, self).__init__(subj, obj, stmt, citation,
                                                     evidence, annotations)
        self.monomer_name = monomer_name
        self.mod = mod
        self.mod_pos = mod_pos
        self.relationship = relationship
        self.activity = activity

    def monomers(self, model):
        monomer = get_create_monomer(model, self.monomer_name)
        site = site_name(self)
        active_state = states[self.mod][1]
        create_site(monomer, site, (active_state,))
        create_site(monomer, self.activity, ('inactive', 'active'))

    def assemble(self, model):
        try:
            kf_activation = model.parameters['kf_activation']
        except KeyError:
            kf_activation = Parameter('kf_activation', 1e-6)
            model.add_component(kf_activation)

        m = model.monomers[self.monomer_name]

        site = site_name(self)
        active_state = states[self.mod][1]

        if self.relationship == 'DirectlyIncreases':
            pre_activity_state = 'inactive'
            post_activity_state = 'active'
        elif self.relationship == 'DirectlyDecreases':
            pre_activity_state = 'active'
            post_activity_state = 'inactive'
        else:
            raise Exception("Invalid modification/activity relationship.")

        r = Rule('%s_%s%s_%s' %
                 (self.monomer_name, site, active_state, self.activity),
                 m(**{site:active_state, self.activity:pre_activity_state}) >>
                 m(**{site:active_state, self.activity:post_activity_state}),
                 kf_activation)
        model.add_component(r)

    def __str__(self):
        return ("ActivityModification(%s, %s, %s, %s, %s)" %
                (self.monomer_name, self.mod, self.mod_pos, self.relationship,
                 self.activity))

class ActivatingSubstitution(Statement):
    def __init__(self, monomer_name, wt_residue, pos, sub_residue, activity,
                 subj, obj, stmt, citation, evidence, annotations):
        super(ActivatingSubstitution, self).__init__(subj, obj, stmt, citation,
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
    def __init__(self, gef_name, gef_activity, ras_name,
                 subj, obj, stmt, citation, evidence, annotations):
        super(RasGef, self).__init__(subj, obj, stmt, citation, evidence,
                                     annotations)
        self.gef_name = gef_name
        self.gef_activity = gef_activity
        self.ras_name = ras_name

    def monomers(self, model):
        gef = get_create_monomer(model, self.gef_name)
        create_site(gef, self.gef_activity, ('inactive','active'))
        ras = get_create_monomer(model, self.ras_name)
        create_site(ras, 'GtpBound', ('inactive', 'active'))

    def assemble(self, model):
        try:
            kf_gef = model.parameters['kf_gef']
        except KeyError:
            kf_gef = Parameter('kf_gef', 1e-6)
            model.add_component(kf_gef)

        gef = model.monomers[self.gef_name]
        ras = model.monomers[self.ras_name]

        r = Rule('%s_activates_%s' %
                 (self.gef_name, self.ras_name),
                 gef(**{self.gef_activity:'active'}) +
                 ras(**{'GtpBound':'inactive'}) >>
                 gef(**{self.gef_activity:'active'}) +
                 ras(**{'GtpBound':'active'}),
                 kf_gef)
        model.add_component(r)

    def __str__(self):
        return ("RasGef(%s, %s, %s)" %
                (self.gef_name, self.gef_activity, self.ras_name))

class RasGap(Statement):
    def __init__(self, gap_name, gap_activity, ras_name,
                 subj, obj, stmt, citation, evidence, annotations):
        super(RasGap, self).__init__(subj, obj, stmt, citation, evidence,
                                     annotations)
        self.gap_name = gap_name
        self.gap_activity = gap_activity
        self.ras_name = ras_name

    def monomers(self, model):
        gap = get_create_monomer(model, self.gap_name)
        create_site(gap, self.gap_activity, ('inactive','active'))
        ras = get_create_monomer(model, self.ras_name)
        create_site(ras, 'GtpBound', ('inactive', 'active'))

    def assemble(self, model):
        try:
            kf_gap = model.parameters['kf_gap']
        except KeyError:
            kf_gap = Parameter('kf_gap', 1e-6)
            model.add_component(kf_gap)

        gap = model.monomers[self.gap_name]
        ras = model.monomers[self.ras_name]

        r = Rule('%s_inactivates_%s' %
                 (self.gap_name, self.ras_name),
                 gap(**{self.gap_activity:'active'}) +
                 ras(**{'GtpBound':'active'}) >>
                 gap(**{self.gap_activity:'active'}) +
                 ras(**{'GtpBound':'inactive'}),
                 kf_gap)
        model.add_component(r)

    def __str__(self):
        return ("RasGap(%s, %s, %s)" %
                (self.gap_name, self.gap_activity, self.ras_name))



