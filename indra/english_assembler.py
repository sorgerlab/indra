import indra.statements as ist

class EnglishAssembler(object):
    def __init__(self, stmts=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        self.model = None

    def add_statements(self, stmts):
        self.statements += stmts

    def make_model(self):
        stmt_strs = []
        for stmt in self.statements:
            if isinstance(stmt, ist.Phosphorylation):
                stmt_strs.append(assemble_phosphorylation(stmt))
            elif isinstance(stmt, ist.Dephosphorylation):
                stmt_strs.append(assemble_dephosphorylation(stmt))
            elif isinstance(stmt, ist.Complex):
                stmt_strs.append(assemble_complex(stmt))
            else:
                print 'Unhandled statement type.'
        model = ' '.join(stmt_strs)
        return model

def assemble_agent_str(agent):
    agent_str = agent.name
    if not agent.mods and not agent.bound_conditions:
        return agent_str

    # Handle bound conditions
    bound_to = [bc.agent.name for bc in
                agent.bound_conditions if bc.is_bound]
    not_bound_to = [bc.agent.name for bc in
                agent.bound_conditions if not bc.is_bound]
    if bound_to:
        agent_str += ' bound to ' + join_list(bound_to)
        if not_bound_to:
            agent_str += ' and not bound to ' +\
                join_list(not_bound_to)
    else:
        if not_bound_to:
            agent_str += ' not bound to ' +\
                join_list(not_bound_to)

    # Handle modification conditions
    # TODO: handle non-phosphorylation mods
    if agent.mods:
        # Special case
        if len(agent.mods) == 1 and agent.mod_sites[0] is None:
            agent_str = abbrev_prefix[agent.mods[0]] + ' ' + agent_str
        else:
            if agent.bound_conditions:
                agent_str += ' and'
            agent_str += ' phosphorylated on '
            mod_lst = []
            for m, p in zip(agent.mods, agent.mod_sites):
                if p is None:
                    if abbrev_word[m]:
                        mod_lst.append(abbrev_word[m])
                else:
                    mod_lst.append(abbrev_letter[m] + str(p))
            agent_str += join_list(mod_lst)

    return agent_str

def join_list(lst):
    if len(lst) > 2:
        s = ', '.join(lst[:-1])
        s += ' and ' + lst[-1]
    elif len(lst) == 2:
        s = lst[0] + ' and ' + lst[1]
    elif len(lst) == 1:
        s = lst[0]
    else:
        s = ''
    return s

def assemble_phosphorylation(stmt):
    sub_str = assemble_agent_str(stmt.sub)
    if stmt.enz is not None:
        enz_str = assemble_agent_str(stmt.enz)
        stmt_str = enz_str + ' phosphorylates ' + sub_str
    else:
        stmt_str = sub_str + ' is phosphorylated'
    # TODO: mod, mod_pos
    if stmt.mod_pos is None:
        if stmt.mod != 'Phosphorylation':
            mod_str = 'on ' + abbrev_word[stmt.mod]
        else:
            mod_str = ''
    else:
        mod_str = 'on ' + abbrev_letter[stmt.mod] + str(stmt.mod_pos)
    stmt_str += ' ' + mod_str
    return make_sentence(stmt_str)

def assemble_dephosphorylation(stmt):
    sub_str = assemble_agent_str(stmt.sub)
    if stmt.enz is not None:
        enz_str = assemble_agent_str(stmt.enz)
        stmt_str = enz_str + ' dephosphorylates ' + sub_str
    else:
        stmt_str = sub_str + ' is dephosphorylated'
    # TODO: mod, mod_pos
    if stmt.mod_pos is None:
        if stmt.mod != 'Phosphorylation':
            mod_str = 'on ' + abbrev_word[stmt.mod]
        else:
            mod_str = ''
    else:
        mod_str = 'on ' + abbrev_letter[stmt.mod] + str(stmt.mod_pos)
    stmt_str += ' ' + mod_str
    return make_sentence(stmt_str)

def assemble_complex(stmt):
    member_strs = [assemble_agent_str(m) for m in stmt.members]
    stmt_str = member_strs[0] + ' binds ' + join_list(member_strs[1:])
    return make_sentence(stmt_str)

def make_sentence(txt):
    #Make sure first letter is capitalized
    txt = txt.strip(' ')
    txt = txt[0].upper() + txt[1:] + '.'
    return txt

abbrev_letter = {
    'Phosphorylation': 'an unknown site',
    'PhosphorylationSerine': 'S',
    'PhosphorylationThreonine': 'T',
    'PhosphorylationTyrosine': 'Y',
}

abbrev_word = {
    'Phosphorylation': 'an unknown site',
    'PhosphorylationSerine': 'a serine residue',
    'PhosphorylationThreonine': 'a threonine residue',
    'PhosphorylationTyrosine': 'a tyrosine residue',
}

abbrev_prefix = {
    'Phosphorylation': 'phosphorylated',
    'PhosphorylationSerine': 'serine-phosphorylated',
    'PhosphorylationThreonine': 'threonine-phosphorylated',
    'PhosphorylationTyrosine': 'tyrosine-phosphorylated',
}
