from collections import defaultdict

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, get_statement_by_name


def get_row_data(stmt_list, ev_totals=None):
    def name(agent):
        return 'None' if agent is None else agent.name

    stmt_rows = {}
    stmt_counts = defaultdict(lambda: 0)
    argument_counts = defaultdict(lambda: 0)
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
        if ev_totals is None:
            stmt_counts[key] += len(s.evidence)
            argument_counts[key[1:]] += len(s.evidence)
        else:
            stmt_counts[key] += ev_totals[s.get_hash()]
            argument_counts[key[1:]] += ev_totals[s.get_hash()]

    # Sort the rows by count and agent names.
    def process(tpl):
        key, stmts = tpl
        verb = key[0]
        inps = key[1:]
        arg_count = argument_counts[inps]
        sub_count = stmt_counts[key]
        new_key = (arg_count, inps)
        new_key += (sub_count, verb)
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

    inp = key[1]
    StmtClass = get_statement_by_name(verb)
    if verb == 'Complex':
        stmt = StmtClass([make_agent(name) for name in inp])
    elif verb == 'Conversion':
        stmt = StmtClass(make_agent(inp[0]),
                         [make_agent(name) for name in inp[1]],
                         [make_agent(name) for name in inp[2]])
    elif verb == 'ActiveForm' or verb == 'HasActivity':
        stmt = StmtClass(make_agent(inp[0]), inp[1], inp[2])
    else:
        stmt = StmtClass(*[make_agent(name) for name in inp])
    ea = EnglishAssembler([stmt])
    return ea.make_model()[:-1]
