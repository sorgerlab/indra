import logging
from collections import defaultdict
from itertools import permutations
from numpy import array, concatenate, zeros

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, Influence, Event, get_statement_by_name


logger = logging.getLogger(__name__)


def _get_keyed_stmts(stmt_list):
    def name(agent):
        return 'None' if agent is None else agent.name

    for s in stmt_list:
        # Create a key.
        verb = s.__class__.__name__
        key = (verb,)
        ags = s.agent_list()
        if verb == 'Complex':
            ag_ns = {name(ag) for ag in ags}
            if 1 < len(ag_ns) < 6:
                for pair in permutations(ag_ns, 2):
                    yield key + tuple(pair),  s
            if len(ag_ns) == 2:
                continue
            key += tuple(sorted(ag_ns))
        elif verb == 'Conversion':
            subj = name(s.subj)
            objs_from = {name(ag) for ag in s.obj_from}
            objs_to = {name(ag) for ag in s.obj_to}
            key += (subj, tuple(sorted(objs_from)), tuple(sorted(objs_to)))
        elif verb == 'ActiveForm':
            key += (name(ags[0]), s.activity, s.is_active)
        elif verb == 'HasActivity':
            key += (name(ags[0]), s.activity, s.has_activity)
        elif verb == 'Influence':
            sns, sid = s.subj.concept.get_grounding()
            ons, oid = s.obj.concept.get_grounding()
            skey = s.subj.concept.name if not sid \
                else sid.split('/')[-1].replace('_', ' ')
            okey = s.obj.concept.name if not oid \
                else oid.split('/')[-1].replace('_', ' ')
            key += (skey, okey)
        else:
            key += tuple([name(ag) for ag in ags])

        yield key, s


def group_and_sort_statements(stmt_list, ev_totals=None, source_counts=None):
    """Group statements by type and arguments, and sort by prevalence.

    Parameters
    ----------
    stmt_list : list[Statement]
        A list of INDRA statements.
    ev_totals : dict{int: int}
        A dictionary, keyed by statement hash (shallow) with counts of total
        evidence as the values. Including this will allow statements to be
        better sorted.
    source_counts : dict{int: OrderedDict}
        A dictionary, keyed by statement hash, with an OrderedDict of
        counts per source (ordered).

    Returns
    -------
    sorted_groups : list[tuple]
        A list of tuples containing a sort key, the statement type, and a list
        of statements, also sorted by evidence count, for that key and type.
        The sort key contains a count of statements with those argument, the
        arguments (normalized strings), the count of statements with those
        arguements and type, and then the statement type.
    """
    if source_counts:
        src_arrays = {}
        for h, counts in source_counts.items():
            src_arrays[h] = array(list(counts.values()))
        n_counts = len(counts) + 1
        cols = list(counts.keys())
    else:
        src_arrays = {}
        n_counts = 1
        cols = None

    def _counts(stmt):
        sh = stmt.get_hash()
        if ev_totals is None or sh not in ev_totals:
            tot = len(stmt.evidence)
        else:
            tot = ev_totals[sh]

        if src_arrays and sh in src_arrays:
            return concatenate([array([tot]), src_arrays[sh]])
        else:
            return array([tot])

    stmt_rows = defaultdict(list)
    stmt_counts = defaultdict(lambda: zeros(n_counts))
    arg_counts = defaultdict(lambda: zeros(n_counts))
    for key, s in _get_keyed_stmts(stmt_list):
        # Update the counts, and add key if needed.
        stmt_rows[key].append(s)

        # Keep track of the total evidence counts for this statement and the
        # arguments.
        stmt_counts[key] += _counts(s)

        # Add up the counts for the arguments, pairwise for Complexes and
        # Conversions. This allows, for example, a complex between MEK, ERK,
        # and something else to lend weight to the interactions between MEK
        # and ERK.
        if key[0] == 'Conversion':
            subj = key[1]
            for obj in key[2] + key[3]:
                arg_counts[(subj, obj)] += _counts(s)
        else:
            arg_counts[key[1:]] += _counts(s)

    # Sort the rows by count and agent names.
    def processed_rows(stmt_rows):
        for key, stmts in stmt_rows.items():
            verb = key[0]
            inps = key[1:]
            sub_count = stmt_counts[key][0]
            arg_count = arg_counts[inps][0]
            if verb == 'Complex' and sub_count == arg_count and len(inps) <= 2:
                if all([len(set(ag.name for ag in s.agent_list())) > 2
                        for s in stmts]):
                    continue
            new_key = (arg_count, inps, sub_count, verb)
            stmts = sorted(stmts,
                           key=(lambda s:
                                _counts(s)[0] + 1 / (1+len(s.agent_list()))),
                           reverse=True)
            if source_counts:
                yield new_key, verb, stmts,\
                      dict(zip(cols, arg_counts[inps][1:].astype(int))),\
                      dict(zip(cols, stmt_counts[key][1:].astype(int)))
            else:
                yield new_key, verb, stmts

    sorted_groups = sorted(processed_rows(stmt_rows),
                           key=lambda tpl: tpl[0], reverse=True)

    return sorted_groups


def make_stmt_from_sort_key(key, verb, agents=None):
    """Make a Statement from the sort key.

    Specifically, the sort key used by `group_and_sort_statements`.
    """
    def make_agent(name):
        if name == 'None' or name is None:
            return None
        return Agent(name)

    StmtClass = get_statement_by_name(verb)
    inps = list(key[1])
    if agents is None:
        agents = []
    if verb == 'Complex':
        agents.extend([make_agent(name) for name in inps])
        stmt = StmtClass(agents[:])
    elif verb == 'Conversion':
        names_from = [make_agent(name) for name in inps[1]]
        names_to = [make_agent(name) for name in inps[2]]
        agents.extend(names_from + names_to)
        stmt = StmtClass(make_agent(inps[0]), names_from, names_to)
    elif verb == 'ActiveForm' or verb == 'HasActivity':
        agents.extend([make_agent(inps[0])])
        stmt = StmtClass(agents[0], inps[1], inps[2])
    elif verb == 'Influence':
        agents.extend([make_agent(inp) for inp in inps[:2]])
        stmt = Influence(*[Event(ag) for ag in agents])
    else:
        agents.extend([make_agent(name) for name in inps])
        stmt = StmtClass(*agents)
    return stmt


def stmt_to_english(stmt):
    """Return an English assembled Statement as a sentence."""
    ea = EnglishAssembler([stmt])
    return ea.make_model()[:-1]


def make_string_from_sort_key(key, verb):
    """Make a Statement string via EnglishAssembler from the sort key.

    Specifically, the sort key used by `group_and_sort_statements`.
    """
    stmt = make_stmt_from_sort_key(key, verb)
    return stmt_to_english(stmt)


def get_simplified_stmts(stmts):
    simple_stmts = []
    for key, s in _get_keyed_stmts(stmts):
        simple_stmts.append(make_stmt_from_sort_key(key, s.__class__.__name__))
    return simple_stmts


def _str_conversion_bits(tpl):
    bolds = ['<b>%s</b>' % el for el in tpl]
    return ', '.join(bolds[:-1]) + ', and ' + bolds[-1]


def make_top_level_label_from_names_key(names):
    try:
        if len(names) == 3 and isinstance(names[1], tuple):  # Conversions
            el_from = _str_conversion_bits(names[1])
            el_to = _str_conversion_bits(names[2])
            tl_label = ("<b>%s</b> converts %s to %s"
                        % (names[0], el_from, el_to))
        else:
            b_names = ['<b>%s</b>' % name for name in names]
            if len(names) == 1:
                tl_label = b_names[0]
            elif len(names) == 2:  # Singleton Modifications
                if names[0] is None or names[0] == 'None':
                    tl_label = b_names[1] + " is modified"
                else:
                    tl_label = b_names[0] + " affects " + b_names[1]
            elif names[1] == "activity":  # ActiveForms
                if names[2] or names[2] == "True":
                    tl_label = b_names[0] + " is active"
                else:
                    tl_label = b_names[0] + " is not active"
            else:  # Large Complexes
                tl_label = b_names[0] + " affects "
                tl_label += ", ".join(b_names[1:-1]) + ', and ' + b_names[-1]
        return tl_label
    except Exception as e:
        logger.error("Could not handle: %s" % str(names))
        raise e
