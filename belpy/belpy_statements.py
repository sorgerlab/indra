
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
        # A default parameter object for phosphorylation
        kf_phospho = Parameter('kf_phospho', 1.)
        model.add_component(kf_phospho)

    def __repr__(self):
        return ("Phosphorylation(%s, %s, %s, %s, %s, %s, %s)" %
                (self.kin_name, self.sub_name, self.mod, self.mod_pos,
                 self.subj, self.obj, self.stmt))

    def __str__(self):
        return ("Phosphorylation(%s, %s, %s, %s)" %
                (self.kin_name, self.sub_name, self.mod, self.mod_pos))

class ActivatingModification(Statement):
    def __init__(self, monomer_name, mod_site, mod_pos, mod_state, activity,
                 subj, obj, stmt):
        super(ActivatingModification, self).__init__(subj, obj, stmt)
        self.monomer_name = monomer_name
        self.mod_site = mod_site
        self.mod_pos = mod_pos
        self.mod_state = mod_state
        self.activity = activity

    def assemble(self, model):
        pass

    def __repr__(self):
        return ("ActivatingModification(%s, %s, %s, %s, %s, %s, %s, %s)" %
                (self.monomer_name, self.mod_site, self.mod_pos,
                 self.mod_state, self.activity, self.subj, self.obj,
                 self.stmt))

    def __str__(self):
        return ("ActivatingModification(%s, %s, %s, %s, %s)" %
                (self.monomer_name, self.mod_site, self.mod_pos,
                 self.mod_state, self.activity))

