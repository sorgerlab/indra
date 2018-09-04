from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import indra.statements as ist

logger = logging.getLogger('english_assembler')

class EnglishAssembler(object):
    """This assembler generates English sentences from INDRA Statements.

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be added to the assembler.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to assemble.
    model : str
        The assembled sentences as a single string.
    """
    def __init__(self, stmts=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        self.model = None

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        self.statements += stmts

    def make_model(self):
        """Assemble text from the set of collected INDRA Statements.

        Returns
        -------
        stmt_strs : str
            Return the assembled text as unicode string. By default, the text
            is a single string consisting of one or more sentences with
            periods at the end.
        """
        stmt_strs = []
        for stmt in self.statements:
            if isinstance(stmt, ist.Modification):
                stmt_strs.append(_assemble_modification(stmt))
            elif isinstance(stmt, ist.Autophosphorylation):
                stmt_strs.append(_assemble_autophosphorylation(stmt))
            elif isinstance(stmt, ist.Complex):
                stmt_strs.append(_assemble_complex(stmt))
            elif isinstance(stmt, ist.RegulateActivity):
                stmt_strs.append(_assemble_regulate_activity(stmt))
            elif isinstance(stmt, ist.RegulateAmount):
                stmt_strs.append(_assemble_regulate_amount(stmt))
            elif isinstance(stmt, ist.ActiveForm):
                stmt_strs.append(_assemble_activeform(stmt))
            elif isinstance(stmt, ist.Translocation):
                stmt_strs.append(_assemble_translocation(stmt))
            elif isinstance(stmt, ist.Gef):
                stmt_strs.append(_assemble_gef(stmt))
            elif isinstance(stmt, ist.Gap):
                stmt_strs.append(_assemble_gap(stmt))
            else:
                logger.warning('Unhandled statement type: %s.' % type(stmt))
        if stmt_strs:
            return ' '.join(stmt_strs)
        else:
            return ''

def _assemble_agent_str(agent):
    """Assemble an Agent object to text."""
    agent_str = agent.name

    # Only do the more detailed assembly for molecular agents
    if not isinstance(agent, ist.Agent):
        return agent_str

    # Handle mutation conditions
    if agent.mutations:
        mut_strs = []
        for mut in agent.mutations:
            res_to = mut.residue_to if mut.residue_to else ''
            res_from = mut.residue_from if mut.residue_from else ''
            pos = mut.position if mut.position else ''
            mut_str = '%s%s%s' % (res_from, pos, res_to)
            mut_strs.append(mut_str)
        mut_strs = '/'.join(mut_strs)
        agent_str = '%s-%s' % (agent_str, mut_strs)

    # Handle location
    if agent.location is not None:
        agent_str += ' in the ' + agent.location

    if not agent.mods and not agent.bound_conditions and not agent.activity:
        return agent_str

    # Handle bound conditions
    bound_to = [bc.agent.name for bc in
                agent.bound_conditions if bc.is_bound]
    not_bound_to = [bc.agent.name for bc in
                agent.bound_conditions if not bc.is_bound]
    if bound_to:
        agent_str += ' bound to ' + _join_list(bound_to)
        if not_bound_to:
            agent_str += ' and not bound to ' +\
                _join_list(not_bound_to)
    else:
        if not_bound_to:
            agent_str += ' not bound to ' +\
                _join_list(not_bound_to)

    # Handle modification conditions
    if agent.mods:
        # Special case
        if len(agent.mods) == 1 and agent.mods[0].position is None:
            prefix = _mod_state_str(agent.mods[0].mod_type)
            if agent.mods[0].residue is not None:
                residue_str =\
                    ist.amino_acids[agent.mods[0].residue]['full_name']
                prefix = residue_str + '-' + prefix
            agent_str =  prefix + ' ' + agent_str
        else:
            if agent.bound_conditions:
                agent_str += ' and'
            agent_str += ' %s on ' % _mod_state_str(agent.mods[0].mod_type)
            mod_lst = []
            for m in agent.mods:
                if m.position is None:
                    if m.residue is not None:
                        residue_str =\
                            ist.amino_acids[m.residue]['full_name']
                        mod_lst.append(residue_str)
                    else:
                        mod_lst.append('an unknown residue')
                elif m.position is not None and m.residue is None:
                    mod_lst.append('amino acid %s' % m.position)
                else:
                    mod_lst.append(m.residue + m.position)
            agent_str += _join_list(mod_lst)


    # Handle activity conditions
    if agent.activity is not None:
        # TODO: handle activity types
        if agent.activity.is_active:
            prefix = 'active'
        else:
            prefix = 'inactive'
        agent_str = prefix + ' ' + agent_str

    return agent_str

def _join_list(lst):
    """Join a list of words in a gramatically correct way."""
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

def _assemble_activeform(stmt):
    """Assemble ActiveForm statements into text."""
    subj_str = _assemble_agent_str(stmt.agent)
    if stmt.is_active:
        is_active_str = 'active'
    else:
        is_active_str = 'inactive'
    if stmt.activity == 'activity':
        stmt_str = subj_str + ' is ' + is_active_str
    elif stmt.activity == 'kinase':
        stmt_str = subj_str + ' is kinase-' + is_active_str
    elif stmt.activity == 'phosphatase':
        stmt_str = subj_str + ' is phosphatase-' + is_active_str
    elif stmt.activity == 'catalytic':
        stmt_str = subj_str + ' is catalytically ' + is_active_str
    elif stmt.activity == 'transcription':
        stmt_str = subj_str + ' is transcriptionally ' + is_active_str
    elif stmt.activity == 'gtpbound':
        stmt_str = subj_str + ' is GTP-bound ' + is_active_str
    return _make_sentence(stmt_str)

def _assemble_modification(stmt):
    """Assemble Modification statements into text."""
    sub_str = _assemble_agent_str(stmt.sub)
    if stmt.enz is not None:
        enz_str = _assemble_agent_str(stmt.enz)
        if _get_is_direct(stmt):
            mod_str = ' ' + _mod_process_verb(stmt) + ' '
        else:
            mod_str = ' leads to the ' + _mod_process_noun(stmt) + ' of '
        stmt_str = enz_str + mod_str + sub_str
    else:
        stmt_str = sub_str + ' is ' + _mod_state_stmt(stmt)

    if stmt.residue is not None:
        if stmt.position is None:
            mod_str = 'on ' + ist.amino_acids[stmt.residue]['full_name']
        else:
            mod_str = 'on ' + stmt.residue + stmt.position
    else:
        mod_str = ''
    stmt_str += ' ' + mod_str
    return _make_sentence(stmt_str)

def _assemble_complex(stmt):
    """Assemble Complex statements into text."""
    member_strs = [_assemble_agent_str(m) for m in stmt.members]
    stmt_str = member_strs[0] + ' binds ' + _join_list(member_strs[1:])
    return _make_sentence(stmt_str)

def _assemble_autophosphorylation(stmt):
    """Assemble Autophosphorylation statements into text."""
    enz_str = _assemble_agent_str(stmt.enz)
    stmt_str = enz_str + ' phosphorylates itself'
    if stmt.residue is not None:
        if stmt.position is None:
            mod_str = 'on ' + ist.amino_acids[stmt.residue]['full_name']
        else:
            mod_str = 'on ' + stmt.residue + stmt.position
    else:
        mod_str = ''
    stmt_str += ' ' + mod_str
    return _make_sentence(stmt_str)

def _assemble_regulate_activity(stmt):
    """Assemble RegulateActivity statements into text."""
    subj_str = _assemble_agent_str(stmt.subj)
    obj_str = _assemble_agent_str(stmt.obj)
    if stmt.is_activation:
        rel_str = ' activates '
    else:
        rel_str = ' inhibits '
    stmt_str = subj_str + rel_str + obj_str
    return _make_sentence(stmt_str)

def _assemble_regulate_amount(stmt):
    """Assemble RegulateAmount statements into text."""
    obj_str = _assemble_agent_str(stmt.obj)
    if stmt.subj is not None:
        subj_str = _assemble_agent_str(stmt.subj)
        if isinstance(stmt, ist.IncreaseAmount):
            rel_str = ' increases the amount of '
        elif isinstance(stmt, ist.DecreaseAmount):
            rel_str = ' decreases the amount of '
        stmt_str = subj_str + rel_str + obj_str
    else:
        if isinstance(stmt, ist.IncreaseAmount):
            stmt_str = obj_str + ' is produced'
        elif isinstance(stmt, ist.DecreaseAmount):
            stmt_str = obj_str + ' is degraded'
    return _make_sentence(stmt_str)

def _assemble_translocation(stmt):
    """Assemble Translocation statements into text."""
    agent_str = _assemble_agent_str(stmt.agent)
    stmt_str = agent_str + ' translocates'
    if stmt.from_location is not None:
        stmt_str += ' from the ' + stmt.from_location
    if stmt.to_location is not None:
        stmt_str += ' to the ' + stmt.to_location
    return _make_sentence(stmt_str)

def _assemble_gap(stmt):
    """Assemble Gap statements into text."""
    subj_str = _assemble_agent_str(stmt.gap)
    obj_str = _assemble_agent_str(stmt.ras)
    stmt_str = subj_str + ' is a GAP for ' + obj_str
    return _make_sentence(stmt_str)

def _assemble_gef(stmt):
    """Assemble Gef statements into text."""
    subj_str = _assemble_agent_str(stmt.gef)
    obj_str = _assemble_agent_str(stmt.ras)
    stmt_str = subj_str + ' is a GEF for ' + obj_str
    return _make_sentence(stmt_str)

def _make_sentence(txt):
    """Make a sentence from a piece of text."""
    #Make sure first letter is capitalized
    txt = txt.strip(' ')
    txt = txt[0].upper() + txt[1:] + '.'
    return txt

def _get_is_direct(stmt):
    '''Returns true if there is evidence that the statement is a direct
    interaction. If any of the evidences associated with the statement
    indicates a direct interatcion then we assume the interaction
    is direct. If there is no evidence for the interaction being indirect
    then we default to direct.'''
    any_indirect = False
    for ev in stmt.evidence:
        if ev.epistemics.get('direct') is True:
            return True
        elif ev.epistemics.get('direct') is False:
            # This guarantees that we have seen at least
            # some evidence that the statement is indirect
            any_indirect = True
    if any_indirect:
        return False
    return True

def _get_is_hypothesis(stmt):
    '''Returns true if there is evidence that the statement is only
    hypothetical. If all of the evidences associated with the statement
    indicate a hypothetical interaction then we assume the interaction
    is hypothetical.'''
    for ev in stmt.evidence:
        if not ev.epistemics.get('hypothesis') is True:
            return True
    return False

def _get_is_hypothesis_adverb(stmt):
    '''Returns the string associated with a statement being hypothetical.'''
    if _get_is_hypothesis(stmt):
        return ' hypothetically '
    else:
        return ''

def _mod_process_verb(stmt):
    mod_name = stmt.__class__.__name__.lower()
    return mod_process_prefix.get(mod_name)

def _mod_process_noun(stmt):
    mod_name = stmt.__class__.__name__.lower()
    return mod_name

def _mod_state_stmt(stmt):
    mod_name = stmt.__class__.__name__.lower()
    return mod_state_prefix.get(mod_name)

def _mod_state_str(s):
    return mod_state_prefix.get(s)

mod_state_prefix = {
    'phosphorylation': 'phosphorylated',
    'dephosphorylation': 'dephosphorylated',
    'ubiquitination': 'ubiquitinated',
    'deubiquitination': 'deubiquitinated',
    'acetylation': 'acetylated',
    'deacetylation': 'deacetylated',
    'hydroxylation': 'hydroxylated',
    'dehydroxylation': 'dehydroxylated',
    'sumoylation': 'sumoylated',
    'desumoylation': 'desumoylated',
    'farnesylation': 'farnesylated',
    'defarnesylation': 'defarnesylated',
    'glycosylation': 'glycosylated',
    'deglycosylation': 'deglycosylated',
    'ribosylation': 'ribosylated',
    'deribosylation': 'deribosylated',
    'modification': 'modified',
}

mod_process_prefix = {
    'phosphorylation': 'phosphorylates',
    'dephosphorylation': 'dephosphorylates',
    'ubiquitination': 'ubiquitinates',
    'deubiquitination': 'deubiquitinates',
    'acetylation': 'acetylates',
    'deacetylation': 'deacetylates',
    'hydroxylation': 'hydroxylates',
    'dehydroxylation': 'dehydroxylates',
    'sumoylation': 'sumoylates',
    'desumoylation': 'desumoylates',
    'farnesylation': 'farnesylates',
    'defarnesylation': 'defarnesylates',
    'glycosylation': 'glycosylates',
    'deglycosylation': 'deglycosylates',
    'ribosylation': 'ribosylates',
    'deribosylation': 'deribosylates'
    }
