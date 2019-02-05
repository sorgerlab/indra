from collections import defaultdict
from itertools import permutations

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, get_statement_by_name


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
            if len(ag_ns) <= 2:
                continue
            key += tuple(sorted(ag_ns))
        elif verb == 'Conversion':
            subj = name(ags[0])
            objs_from = {name(ag) for ag in ags[1]}
            objs_to = {name(ag) for ag in ags[2]}
            key += (subj, tuple(sorted(objs_from)), tuple(sorted(objs_to)))
        elif verb == 'ActiveForm':
            key += (name(ags[0]), s.activity, s.is_active)
        elif verb == 'HasActivity':
            key += (name(ags[0]), s.activity, s.has_activity)
        else:
            key += tuple([name(ag) for ag in ags])

        yield key, s


def get_row_data(stmt_list, ev_totals=None):
    def _count(stmt):
        if ev_totals is None:
            return len(stmt.evidence)
        else:
            return ev_totals[stmt.get_hash()]

    stmt_rows = defaultdict(list)
    stmt_counts = defaultdict(lambda: 0)
    arg_counts = defaultdict(lambda: 0)
    for key, s in _get_keyed_stmts(stmt_list):
        # Update the counts, and add key if needed.
        stmt_rows[key].append(s)

        # Keep track of the total evidence counts for this statement and the
        # arguments.
        stmt_counts[key] += _count(s)

        # Add up the counts for the arguments, pairwise for Complexes and
        # Conversions. This allows, for example, a complex between MEK, ERK,
        # and something else to lend weight to the interactions between MEK
        # and ERK.
        if key[0] == 'Conversion':
            subj = key[1]
            for obj in key[2] + key[3]:
                arg_counts[(subj, obj)] += _count(s)
        else:
            arg_counts[key[1:]] += _count(s)

    # Sort the rows by count and agent names.
    def process_rows(stmt_rows):
        for key, stmts in stmt_rows.items():
            verb = key[0]
            inps = key[1:]
            sub_count = stmt_counts[key]
            arg_count = arg_counts[inps]
            if verb == 'Complex' and sub_count == arg_count and len(inps) == 2:
                continue
            new_key = (arg_count, inps, sub_count, verb)
            yield new_key, verb, sorted(stmts, key=_count, reverse=True)

    row_data = sorted(process_rows(stmt_rows),
                      key=lambda tpl: tpl[0], reverse=True)

    return row_data


def make_statement_string(key, verb):
    """Make a Statement string via EnglishAssembler from `get_row_data` key."""
    def make_agent(name):
        if name == 'None' or name is None:
            return None
        return Agent(name)

    StmtClass = get_statement_by_name(verb)
    inps = list(key[1])
    if verb == 'Complex':
        stmt = StmtClass([make_agent(name) for name in inps])
    elif verb == 'Conversion':
        stmt = StmtClass(make_agent(inps[0]),
                         [make_agent(name) for name in inps[1]],
                         [make_agent(name) for name in inps[2]])
    elif verb == 'ActiveForm' or verb == 'HasActivity':
        stmt = StmtClass(make_agent(inps[0]), inps[1], inps[2])
    else:
        stmt = StmtClass(*[make_agent(name) for name in inps])
    ea = EnglishAssembler([stmt])
    return ea.make_model()[:-1]
