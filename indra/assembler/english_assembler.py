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
    if agent.mods:
        # Special case
        if len(agent.mods) == 1 and agent.mods[0].position is None:
            prefix = abbrev_prefix[agent.mods[0].mod_type]
            if agent.mods[0].residue is not None:
                residue_str =\
                    ist.amino_acids[agent.mods[0].residue]['full_name']
                prefix = residue_str + '-' + prefix
            agent_str =  prefix + ' ' + agent_str
        else:
            if agent.bound_conditions:
                agent_str += ' and'
            agent_str += ' %s on ' % abbrev_prefix[agent.mods[0].mod_type]
            mod_lst = []
            for m in agent.mods:
                if m.position is None:
                    if m.residue is not None:
                        residue_str =\
                            ist.amino_acids[m.residue]['full_name']
                        mod_lst.append(residue_str)
                    else:
                        mod_lst.append('an unknown residue')
                else:
                    mod_lst.append(m.residue + m.position)
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

    if stmt.residue is not None:
        if stmt.position is None:
            mod_str = 'on ' + ist.amino_acids[stmt.residue]['full_name']
        else:
            mod_str = 'on ' + stmt.residue + stmt.position
    else:
        mod_str = ''
    stmt_str += ' ' + mod_str
    return make_sentence(stmt_str)

def assemble_dephosphorylation(stmt):
    sub_str = assemble_agent_str(stmt.sub)
    if stmt.enz is not None:
        enz_str = assemble_agent_str(stmt.enz)
        stmt_str = enz_str + ' dephosphorylates ' + sub_str
    else:
        stmt_str = sub_str + ' is dephosphorylated'

    if stmt.residue is not None:
        if stmt.position is None:
            mod_str = 'on ' + ist.amino_acids[stmt.residue]['full_name']
        else:
            mod_str = 'on ' + stmt.residue + stmt.position
    else:
        mod_str = 'on an unknown residue '
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

abbrev_prefix = {
    'phosphorylation': 'phosphorylated',
    'ubiquitination': 'ubiquitinated'
}
