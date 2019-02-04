from collections import defaultdict
from itertools import permutations

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, get_statement_by_name


def get_row_data(stmt_list, ev_totals=None):
    def name(agent):
        return 'None' if agent is None else agent.name

    stmt_rows = {}
    stmt_counts = defaultdict(lambda: 0)
    arg_counts = defaultdict(lambda: 0)
    for s in stmt_list:
        # Create a key.
        verb = s.__class__.__name__
        key = (verb,)

        ags = s.agent_list()
        if verb == 'Complex':
            ag_ns = {name(ag) for ag in ags}
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

        # Update the counts, and add key if needed.
        if key not in stmt_rows.keys():
            stmt_rows[key] = []
        stmt_rows[key].append(s)

        # Keep track of the total evidence counts for this statement and the
        # arguments.
        if ev_totals is None:
            stmt_counts[key] += len(s.evidence)
            arg_count = len(s.evidence)
        else:
            stmt_counts[key] += ev_totals[s.get_hash()]
            arg_count = ev_totals[s.get_hash()]

        # Add up the counts for the arguments, pairwise for Complexes and
        # Conversions. This allows, for example, a complex between MEK, ERK,
        # and something else to lend weight to the interactions between MEK
        # and ERK.
        if verb == 'Complex' and len(key[1:]) > 1:
            for pair in permutations(key[1:], 2):
                arg_counts[pair] += arg_count
        elif verb == 'Conversion':
            subj = key[1]
            for obj in key[2] + key[3]:
                arg_counts[(subj, obj)] += arg_count
        else:
            arg_counts[key[1:]] += arg_count

    # Sort the rows by count and agent names.
    def process(tpl):
        key, stmts = tpl
        verb = key[0]
        inps = key[1:]
        sub_count = stmt_counts[key]
        if verb == 'Complex' and len(inps) > 1:
            arg_count, key_inps = max((arg_counts[p], p)
                                      for p in permutations(inps, 2))
            new_key = (arg_count, key_inps, sub_count, verb, inps)
        elif verb == 'Conversion':
            subj = inps[0]
            arg_count, key_inps = max((arg_counts[(subj, obj)], (subj, obj))
                                      for obj in inps[1] + inps[2])
            new_key = (arg_count, key_inps, sub_count, verb, inps)
        else:
            arg_count = arg_counts[inps]
            new_key = (arg_count, inps, sub_count, verb)
        return new_key, verb, stmts

    row_data = sorted((process(t) for t in stmt_rows.items()),
                      key=lambda tpl: tpl[0], reverse=True)

    return row_data


def make_statement_string(key, verb):
    """Make a Statement string via EnglishAssembler from `get_row_data` key."""
    def make_agent(name):
        if name == 'None' or name is None:
            return None
        return Agent(name)

    StmtClass = get_statement_by_name(verb)
    if verb == 'Complex':
        inps = key[-1]
        stmt = StmtClass([make_agent(name) for name in inps])
    elif verb == 'Conversion':
        inps = key[-1]
        stmt = StmtClass(make_agent(inps[0]),
                         [make_agent(name) for name in inps[1]],
                         [make_agent(name) for name in inps[2]])
    elif verb == 'ActiveForm' or verb == 'HasActivity':
        inps = key[1]
        stmt = StmtClass(make_agent(inps[0]), inps[1], inps[2])
    else:
        inps = key[1]
        stmt = StmtClass(*[make_agent(name) for name in inps])
    ea = EnglishAssembler([stmt])
    return ea.make_model()[:-1]
