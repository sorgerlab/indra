import logging
from collections import defaultdict
from itertools import permutations
from numpy import array, concatenate, zeros

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, Influence, Event, get_statement_by_name


logger = logging.getLogger(__name__)


db_sources = ['phosphosite', 'cbn', 'pc11', 'biopax', 'bel_lc',
              'signor', 'biogrid', 'lincs_drug', 'tas', 'hprd', 'trrust',
              'ctd', 'virhostnet', 'phosphoelm', 'drugbank', 'omnipath']

reader_sources = ['geneways', 'tees', 'isi', 'trips', 'rlimsp', 'medscan',
                  'sparser', 'eidos', 'reach']


# These are mappings where the actual INDRA source, as it appears
# in the evidence source_api is inconsistent with the colors here and
# with what comes out of the INDRA DB
internal_source_mappings = {
    'bel': 'bel_lc'
}


all_sources = db_sources + reader_sources


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


class StmtMetricAgg:
    """A class to aggregate statement metrics efficiently."""
    def __init__(self, ev_counts=None, beliefs=None,
                 source_count_arrays=None, source_count_cols=None):
        self.__src_arrays = source_count_arrays
        self.__source_cols = source_count_cols
        n_sources = len(next(iter(self.__src_arrays.values())))
        self.__ev_counts = ev_counts
        self.__beliefs = beliefs

        # Init the results.
        self.src_count = zeros(n_sources)
        self.ev_count = 0
        self.belief = 0

    def add(self, stmt):
        sh = stmt.get_hash()

        # Update the evidence count.
        if self.__ev_counts is None or sh not in self.__ev_counts:
            n_ev = len(stmt.evidence)
        else:
            n_ev = self.__ev_counts[sh]
        self.ev_count += n_ev

        # Update the source counts.
        if self.__src_arrays and sh in self.__src_arrays:
            self.src_count += self.__src_arrays[sh]

        # Update the belief.
        if self.__beliefs is None or sh not in self.__beliefs:
            belief = stmt.belief
        else:
            belief = self.__beliefs[sh]
        self.belief = max(self.belief, belief)
        return

    def get_source_counts(self):
        if self.__source_cols is None:
            assert not len(self.src_count)
            return {}
        return dict(zip(self.__source_cols, self.src_count.astype(int)))

    @classmethod
    def get_prepped_constructor(cls, ev_counts=None, beliefs=None,
                                source_counts=None):
        # Turn the source dicts into arrays.
        if source_counts:
            src_arrays = {}
            for h, counts in source_counts.items():
                src_arrays[h] = array(list(counts.values()))
            source_cols = list(counts.keys())
        else:
            src_arrays = {}
            source_cols = None

        # Define the constructor function.
        def constructor():
            return cls(ev_counts, beliefs, src_arrays, source_cols)
        return constructor


def group_and_sort_statements(stmt_list, sort_by='ev_count', ev_counts=None,
                              source_counts=None, beliefs=None):
    """Group statements by type and arguments, and sort by prevalence.

    Parameters
    ----------
    stmt_list : list[Statement]
        A list of INDRA statements.
    sort_by : str
        Indicate which parameter to sort by, either 'belief' or 'ev_count'. The
        default is 'ev_count'.
    ev_counts : dict{int: int}
        A dictionary, keyed by statement hash (shallow) with counts of total
        evidence as the values. Including this will allow statements to be
        sorted by evidence count..
    beliefs : dict{int: float}
        A dictionary, keyed by statement hash (shallow) with the belief of each
        statement. Including this will allow statements to be sorted by belief.
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
        arguments and type, and then the statement type.
    """
    if sort_by not in ['belief', 'ev_count']:
        raise ValueError(f"Parameter `sort_by` must be 'belief' or 'ev_count'!")

    stmt_rows = defaultdict(list)

    counter_maker = StmtMetricAgg\
        .get_prepped_constructor(ev_counts, beliefs, source_counts)
    stmt_counts = defaultdict(counter_maker)
    arg_counts = defaultdict(counter_maker)
    for key, s in _get_keyed_stmts(stmt_list):
        # Update the counts, and add key if needed.
        stmt_rows[key].append(s)

        # Keep track of the evidence counts, source counts, and max belief for
        # this statement and the arguments.
        stmt_counts[key].add(s)

        # Add up the counts for the arguments, pairwise for Complexes and
        # Conversions. This allows, for example, a complex between MEK, ERK,
        # and something else to lend weight to the interactions between MEK
        # and ERK.
        if key[0] == 'Conversion':
            subj = key[1]
            for obj in key[2] + key[3]:
                arg_counts[(subj, obj)].add(s)
        else:
            arg_counts[key[1:]].add(s)

    def _sort_key_func(stmt):
        sh = stmt.get_hash()
        n_agents = len(s.agent_list())
        if sort_by == 'ev_count':
            if ev_counts and sh in ev_counts:
                cnt = ev_counts[sh]
            else:
                cnt = len(stmt.evidence)
            return cnt + 1/(1 + n_agents)
        elif sort_by == 'belief':
            if beliefs and sh in beliefs:
                b = beliefs[sh]
            else:
                b = stmt.belief
            return b + 0.01/(1 + n_agents)
        else:
            assert False, "Invalid sort_by, should not be possible here."

    # Sort the rows by count and agent names.
    def processed_rows(stmt_rows):
        for key, stmts in stmt_rows.items():
            verb = key[0]
            inps = key[1:]
            sub_count = stmt_counts[key].ev_count
            arg_count = arg_counts[inps].ev_count
            if verb == 'Complex' and sub_count == arg_count and len(inps) <= 2:
                if all([len(set(ag.name for ag in s.agent_list())) > 2
                        for s in stmts]):
                    continue
            new_key = (arg_count, inps, sub_count, verb)
            stmts = sorted(stmts, key=_sort_key_func, reverse=True)
            if source_counts:
                yield new_key, verb, stmts, \
                      arg_counts[inps].get_source_counts(), \
                      stmt_counts[key].get_source_counts()
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
    elif verb == 'Association':
        agents.extend([make_agent(inp) for inp in inps])
        stmt = StmtClass([Event(ag) for ag in agents])
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


def standardize_counts(counts):
    """Standardize hash-based counts dicts to be int-keyed."""
    standardized_counts = {}
    for k, v in counts.items():
        try:
            int_k = int(k)
            standardized_counts[int_k] = v
        except ValueError:
            logger.warning('Could not convert statement hash %s to int' % k)
    return standardized_counts


def get_available_ev_counts(stmts):
    return {stmt.get_hash(): len(stmt.evidence) for stmt in stmts}


def get_available_beliefs(stmts):
    return {stmt.get_hash(): stmt.belief for stmt in stmts}


def get_available_source_counts(stmts):
    return {stmt.get_hash(): _get_available_ev_source_counts(stmt.evidence)
            for stmt in stmts}


def _get_available_ev_source_counts(evidences):
    counts = _get_initial_source_counts()
    for ev in evidences:
        sa = internal_source_mappings.get(ev.source_api, ev.source_api)
        try:
            counts[sa] += 1
        except KeyError:
            continue
    return counts


def _get_initial_source_counts():
    return {s: 0 for s in all_sources}
