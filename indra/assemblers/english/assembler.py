from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import indra.statements as ist
from indra.explanation.reporting import PybelEdge

logger = logging.getLogger(__name__)


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
            elif isinstance(stmt, ist.Association):
                stmt_strs.append(_assemble_association(stmt))
            elif isinstance(stmt, ist.Complex):
                stmt_strs.append(_assemble_complex(stmt))
            elif isinstance(stmt, ist.Influence):
                stmt_strs.append(_assemble_influence(stmt))
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
            elif isinstance(stmt, ist.Conversion):
                stmt_strs.append(_assemble_conversion(stmt))
            elif isinstance(stmt, PybelEdge):
                stmt_strs.append(_assemble_pybel_edge(stmt))
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
        is_generic = False
        mut_strs = []
        for mut in agent.mutations:
            res_to = mut.residue_to if mut.residue_to else ''
            res_from = mut.residue_from if mut.residue_from else ''
            pos = mut.position if mut.position else ''
            mut_str = '%s%s%s' % (res_from, pos, res_to)
            # If this is the only mutation and there are no details
            # then this is a generic mutant
            if not mut_str and len(agent.mutations) == 1:
                is_generic = True
                break
            mut_strs.append(mut_str)
        if is_generic:
            agent_str = 'mutated ' + agent_str
        else:
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
            agent_str = prefix + ' ' + agent_str
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
        # Get the modifier specific to the activity type, if any
        pre_prefix = \
            activity_type_prefix.get(agent.activity.activity_type, '')
        if agent.activity.is_active:
            prefix = pre_prefix + 'active'
        else:
            # See if there is a special override for the inactive form
            if agent.activity.activity_type in inactivity_type_prefix_override:
                pre_prefix = inactivity_type_prefix_override[
                    agent.activity.activity_type]
            prefix = pre_prefix + 'inactive'
        agent_str = prefix + ' ' + agent_str

    return agent_str


def english_join(lst):
    """Join a list of strings according to English grammar.

    Parameters
    ----------
    lst : list of str
        A list of strings to join.

    Returns
    -------
    str
        A string which describes the list of elements, e.g.,
        "apples, pears, and bananas".
    """
    return _join_list(lst, oxford=True)


def _join_list(lst, oxford=True):
    """Join a list of words in a gramatically correct way."""
    if len(lst) > 2:
        s = ', '.join(lst[:-1])
        if oxford:
            s += ','
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
    elif stmt.position is not None:
        mod_str = 'at position %s' % stmt.position
    else:
        mod_str = ''
    stmt_str += ' ' + mod_str
    return _make_sentence(stmt_str)


def _assemble_association(stmt):
    """Assemble Association statements into text."""
    member_strs = [_assemble_agent_str(m.concept) for m in stmt.members]
    stmt_str = member_strs[0] + ' is associated with ' + \
        _join_list(member_strs[1:])
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


def _assemble_conversion(stmt):
    """Assemble a Conversion statement into text."""
    reactants = _join_list([_assemble_agent_str(r) for r in stmt.obj_from])
    products = _join_list([_assemble_agent_str(r) for r in stmt.obj_to])

    if stmt.subj is not None:
        subj_str = _assemble_agent_str(stmt.subj)
        stmt_str = '%s catalyzes the conversion of %s into %s' % \
            (subj_str, reactants, products)
    else:
        stmt_str = '%s is converted into %s' % (reactants, products)
    return _make_sentence(stmt_str)


def _assemble_influence(stmt):
    """Assemble an Influence statement into text."""
    subj_str = _assemble_agent_str(stmt.subj.concept)
    obj_str = _assemble_agent_str(stmt.obj.concept)

    # Note that n is prepended to increase to make it "an increase"
    if stmt.subj.delta.polarity is not None:
        subj_delta_str = ' decrease' if stmt.subj.delta.polarity == -1 \
            else 'n increase'
        subj_str = 'a%s in %s' % (subj_delta_str, subj_str)

    if stmt.obj.delta.polarity is not None:
        obj_delta_str = ' decrease' if stmt.obj.delta.polarity == -1 \
            else 'n increase'
        obj_str = 'a%s in %s' % (obj_delta_str, obj_str)

    stmt_str = '%s causes %s' % (subj_str, obj_str)
    return _make_sentence(stmt_str)


def _assemble_pybel_edge(pybel_edge):
    source_str = _assemble_agent_str(pybel_edge.source)
    target_str = _assemble_agent_str(pybel_edge.target)
    if pybel_edge.relation == 'hasComponent':
        if pybel_edge.reverse:
            rel_str = ' is a part of '
        else:
            rel_str = ' has a component '
    elif pybel_edge.relation == 'hasVariant':
        if pybel_edge.reverse:
            rel_str = ' is a variant of '
        else:
            rel_str = ' has a variant '
    edge_str = source_str + rel_str + target_str
    return _make_sentence(edge_str)


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
    """Returns true if there is evidence that the statement is only
    hypothetical. If all of the evidences associated with the statement
    indicate a hypothetical interaction then we assume the interaction
    is hypothetical."""
    for ev in stmt.evidence:
        if not ev.epistemics.get('hypothesis') is True:
            return True
    return False


def _get_is_hypothesis_adverb(stmt):
    """Return the string associated with a statement being hypothetical."""
    if _get_is_hypothesis(stmt):
        return ' hypothetically '
    else:
        return ''


def _mod_process_verb(stmt):
    # Example: Phosphorylation -> phosphorylates
    mod_name = stmt.__class__.__name__.lower()
    return statement_present_verb(mod_name)


def _mod_process_noun(stmt):
    # Example: Phosphorylation -> phosphorylation
    mod_name = stmt.__class__.__name__.lower()
    return mod_name


def _mod_state_stmt(stmt):
    # Example: Phosphorylation -> phosphorylated
    mod_name = stmt.__class__.__name__.lower()
    return statement_passive_verb(mod_name)


def _mod_state_str(s):
    return statement_passive_verb(s)


def statement_passive_verb(stmt_type):
    """Return the passive / state verb form of a statement type.

    Parameters
    ----------
    stmt_type : str
        The lower case string form of a statement type, for instance,
        'phosphorylation'.

    Returns
    -------
    str
        The passive/state verb form of a statement type, for instance,
        'phosphorylated'.
    """
    override = {
        'complex': 'bound',
        'regulateamount': 'amount regulated',
        'decreaseamount': 'decreased',
        'increaseamount': 'increased',
        'gap': 'GAP-regulated',
        'gef': 'GEF-regulated',
        'gtpactivation': 'GTP-activated',
        'influence': 'influenced',
        'event': 'happened',
        'conversion': 'converted',
        'modification': 'modified',
        'addmodification': 'modified',
        'removemodification': 'unmodified',
        'regulateactivity': 'activity regulated',
    }
    return override.get(stmt_type) if stmt_type in override else \
        stmt_type[:-3] + 'ed'


def statement_present_verb(stmt_type):
    """Return the present verb form of a statement type.

    Parameters
    ----------
    stmt_type : str
        The lower case string form of a statement type, for instance,
        'phosphorylation'.

    Returns
    -------
    str
        The present verb form of a statement type, for instance,
        'phosphorylates'.
    """
    override = {
        'complex': 'binds',
        'regulateamount': 'regulates the amount of',
        'increaseamount': 'increases the amount of',
        'decreaseamount': 'decreases the amount of',
        'gef': 'acts as a GEF for',
        'gap': 'acts as a GAP for',
        'inhibition': 'inhibits',
        'gtpactivation': 'activates when bound to GTP',
        'regulateactivity': 'regulates the activity of',
        'activeform': 'has active form',
        'conversion': 'converts',
        'influence': 'influences',
        'modification': 'modifies',
        'addmodification': 'adds a modification to',
        'removemodification': 'removes a modification of',
        'selfmodification': 'modifies itself',
        'event': 'happens'
    }
    return override.get(stmt_type) if stmt_type in override else \
        stmt_type[:-3] + 'es'


def statement_base_verb(stmt_type):
    """Return the base verb form of a statement type.

    Parameters
    ----------
    stmt_type : str
        The lower case string form of a statement type, for instance,
        'phosphorylation'.

    Returns
    -------
    str
        The base verb form of a statement type, for instance, 'phosphorylate'.
    """
    override = {
        'complex': 'bind',
        'regulateamount': 'regulate the amount of',
        'increaseamount': 'increase the amount of',
        'decreaseamount': 'decrease the amount of',
        'gef': 'act as a GEF for',
        'gap': 'act as a GAP for',
        'inhibition': 'inhibit',
        'gtpactivation': 'activate when bound to GTP',
        'regulateactivity': 'regulate the activity of',
        'activeform': 'have active form',
        'conversion': 'convert',
        'influence': 'influence',
        'modification': 'modify',
        'addmodification': 'add a modification to',
        'removemodification': 'remove a modification of',
        'selfmodification': 'modify itself',
        'event': 'happen'
    }
    return override.get(stmt_type) if stmt_type in override \
        else stmt_type[:-3] + 'e'


activity_type_prefix = {
    'catalytic': 'catalytically ',
    'gap': 'GAP-',
    'gef': 'GEF-',
    'gtpbound': 'GTP-bound ',
    'kinase': 'kinase-',
    'phosphatase': 'phosphatase-',
    'transcription': 'transcriptionally ',
    }


inactivity_type_prefix_override = {
    'gtpbound': 'GDP-bound ',
    }
