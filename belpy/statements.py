
class Statement(object):
    def __init__(self, subj, obj, stmt):
        self.subj = subj
        self.obj = obj
        self.stmt = stmt

class Phosphorylation(Statement):
    def __init__(self, kin_name, sub_name, mod, mod_pos, subj, obj, stmt):
        super(Phosphorylation, self).__init__(subj, obj, stmt)
        self.kin_name = kin_name
        self.sub_name = sub_name
        self.mod = mod
        self.mod_pos = mod_pos

    def assemble(self, model):
        kf_phospho = Parameter('kf_phospho', 1.)
        model.add_component(kf_phospho)

    def __str__(self):
        return ("Phosphorylation(%s, %s, %s, %s)" %
                (self.kin_name, self.sub_name, self.mod, self.mod_pos))

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
        kf_activation = Parameter('kf_activation', 1e5)
        model.add_component(kf_activation)

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
        kf_gef = Parameter('kf_gef', 1.)
        model.add_component(kf_gef)

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



