from pysb import *
from rdf_to_pysb import abbrevs, states

class Statement(object):
    def __init__(self, subj, obj, stmt):
        self.subj = subj
        self.obj = obj
        self.stmt = stmt

class Modification(Statement):
    def __init__(self, enz_name, sub_name, mod, mod_pos, subj, obj, stmt):
        super(Modification, self).__init__(subj, obj, stmt)
        self.enz_name = enz_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def __str__(self):
        return ("%s(%s, %s, %s, %s)" %
                (type(self).__name__, self.enz_name, self.sub_name, self.mod,
                 self.mod_pos))

class Phosphorylation(Modification):
    def assemble(self, model):
        try:
            kf_phospho = model.parameters['kf_phospho']
        except KeyError:
            kf_phospho = Parameter('kf_phospho', 1.)
            model.add_component(kf_phospho)

        enz = model.monomers[self.enz_name]
        sub = model.monomers[self.sub_name]

        site_name = '%s%s' % (abbrevs[self.mod], self.mod_pos)
        r = Rule('%s_phospho_%s_%s' %
                 (self.enz_name, self.sub_name, site_name),
                 enz() + sub(**{site_name:'u'}) >>
                 enz() + sub(**{site_name:'p'}),
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

class Dephosphorylation(Statement):
    def __init__(self, phos_name, sub_name, mod, mod_pos, subj, obj, stmt):
        super(Dephosphorylation, self).__init__(subj, obj, stmt)
        self.phos_name = phos_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def assemble(self, model):
        kf_dephospho = Parameter('kf_dephospho', 1.)
        model.add_component(kf_dephospho)

    def __str__(self):
        return ("Dehosphorylation(%s, %s, %s, %s)" %
                (self.phos_name, self.sub_name, self.mod, self.mod_pos))

class ActivatingModification(Statement):
    def __init__(self, monomer_name, mod, mod_pos, activity,
                 subj, obj, stmt):
        super(ActivatingModification, self).__init__(subj, obj, stmt)
        self.monomer_name = monomer_name
        self.mod = mod
        self.mod_pos = mod_pos
        self.activity = activity

    def assemble(self, model):
        try:
            kf_activation = model.parameters['kf_activation']
        except KeyError:
            kf_activation = Parameter('kf_activation', 1e5)
            model.add_component(kf_activation)

        m = model.monomers[self.monomer_name]

        if self.mod_pos is not None:
            site_name = '%s%s' % (abbrevs[self.mod], self.mod_pos)
        else:
            site_name = '%s' % abbrevs[self.mod]
        active_state = states[self.mod][1]

        r = Rule('%s_%s%s_%s' %
                 (self.monomer_name, site_name, active_state, self.activity),
                 m(**{site_name:active_state, self.activity:'inactive'}) >>
                 m(**{site_name:active_state, self.activity:'active'}),
                 kf_activation)
        model.add_component(r)

    def __str__(self):
        return ("ActivatingModification(%s, %s, %s, %s)" %
                (self.monomer_name, self.mod, self.mod_pos, self.activity))

class ActivatingSubstitution(Statement):
    def __init__(self, monomer_name, wt_residue, pos, sub_residue, activity,
                 subj, obj, stmt):
        super(ActivatingSubstitution, self).__init__(subj, obj, stmt)
        self.monomer_name = monomer_name
        self.wt_residue = wt_residue
        self.pos = pos
        self.sub_residue = sub_residue
        self.activity = activity

    def assemble(self, model):
        kf_activation = Parameter('kf_activation', 1e5)
        model.add_component(kf_activation)

    def __str__(self):
        return ("ActivatingSubstitution(%s, %s, %s, %s, %s)" %
                (self.monomer_name, self.wt_residue, self.pos,
                 self.sub_residue, self.activity))

class RasGef(Statement):
    def __init__(self, gef_name, gef_activity, ras_name,
                 subj, obj, stmt):
        super(RasGef, self).__init__(subj, obj, stmt)
        self.gef_name = gef_name
        self.gef_activity = gef_activity
        self.ras_name = ras_name

    def assemble(self, model):
        try:
            kf_gef = model.parameters['kf_gef']
        except KeyError:
            kf_gef = Parameter('kf_gef', 1.)
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
                 subj, obj, stmt):
        super(RasGap, self).__init__(subj, obj, stmt)
        self.gap_name = gap_name
        self.gap_activity = gap_activity
        self.ras_name = ras_name

    def assemble(self, model):
        kf_gap = Parameter('kf_gap', 1.)
        model.add_component(kf_gap)

    def __str__(self):
        return ("RasGap(%s, %s, %s)" %
                (self.gap_name, self.gap_activity, self.ras_name))



