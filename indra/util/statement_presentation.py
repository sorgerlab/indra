import logging
from collections import defaultdict
from itertools import permutations
from numpy import array, concatenate, zeros

from indra.assemblers.english import EnglishAssembler
from indra.statements import Agent, Influence, Event, get_statement_by_name, \
    Statement

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


def _get_relation_keyed_stmts(stmt_list):
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
                    yield key + tuple(pair), s
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


def merge_to_metric_dict(**kwargs):
    """Merge metric dictionaries into one.

    Merge dictionaries keyed by hash, each of which measures one or more
    metrics of a statement, into a single dict keyed by hash.

    The dictionaries MUST have the same keys (i.e. the same hashes). The values
    of a metric may be itself a dict, in which case the key sets will be
    merged to a flat dict (in this way you can extend a metric dict).
    """
    # If there are no kwargs, just return an empty dict.
    if not kwargs:
        return {}

    # Compile the metric dict.
    hash_set = None
    metric_dict = {}
    for metric_name, metrics in kwargs.items():
        # Check if metrics is None. If so, just skip. This happens if e.g. a
        # default argument is passed along
        if not metrics:
            continue

        # If we don't have a hash set yet, make it. Otherwise assert it matches.
        if hash_set is None:
            hash_set = {int(h) for h in metrics.keys()}
            metric_dict = {h: {} for h in hash_set}
        else:
            assert {int(h) for h in metrics.keys()} == hash_set,\
                "Dictionary key sets do not match."

        # Fill up the metric dict.
        values = None
        for h, val in metrics.items():
            # Check that if this is a dictionary, the values match up.
            if isinstance(val, dict):
                if values is None:
                    values = tuple(val.keys())
                    for k in values:
                        assert k not in metric_dict[h], \
                            f"Metric label {k} from dictionary metric " \
                            f"already in use."
                else:
                    assert tuple(val.keys()) == values, "Values to not equal."

                metric_dict[int(h)].update(val)
            else:
                if values is None:
                    values = (metric_name,)
                    assert metric_name not in metric_dict[h], \
                        f"Metric label {metric_name} from keyword already " \
                        f"in use."
                metric_dict[int(h)][metric_name] = val
    return metric_dict


class Metriker:
    def __init__(self, keys, stmt_metrics, original_types):
        self.__metriks = {}
        self.__keys = keys
        self.__stmt_metrics = stmt_metrics
        self.__original_types = original_types
        self.__filling = True

    def __getitem__(self, item):
        if item not in self.__metriks:
            if self.__filling:
                # Remember, this is passing REFERENCES to the stmt_metrics dict.
                self.__metriks[item] = Metrik(self.__keys, self.__stmt_metrics,
                                              self.__original_types)
            else:
                raise KeyError(f"Key \"{item}\" not found! "
                               f"{self.__class__.__name__} is frozen.")
        return self.__metriks[item]

    def _set_metriks_to_stmt_metrics(self):
        for key, arr in self.__stmt_metrics.items():
            self.__metriks[key] = \
                Metrik.from_array(self.__keys, arr, self.__original_types)
        return

    def freeze(self):
        self.__filling = False
        for metrik in self.__metriks.values():
            metrik.freeze()
        return

    @classmethod
    def from_stmt_list(cls, stmt_list, metric_dict=None):
        keys = ('ev_count', 'belief', 'ag_count')
        original_types = (int, float, int)
        if metric_dict:
            if len(stmt_list) != len(metric_dict):
                raise ValueError("The `stmt_list` and `metric_dict` must be "
                                 "the same length.")
            first_metric = next(iter(metric_dict.values()))
            metric_dict_keys = set(first_metric.keys())
            stmt_loop = any(k not in metric_dict_keys for k in keys)
            original_types += tuple(type(first_metric[k])
                                    for k in metric_dict_keys if k not in keys)
            keys += tuple(k for k in metric_dict_keys if k not in keys)
        else:
            stmt_loop = True

        # Populate the dict of arrays, stmt_metrics
        stmt_metrics = {}
        if stmt_loop:

            # Iterate over statements, filling in values that may have been
            # missing.
            for stmt in stmt_list:
                sh = stmt.get_hash()
                values = []
                for k in keys:
                    if metric_dict and k in metric_dict[sh]:
                        values.append(metric_dict[sh][k])
                    elif k == 'ev_count':
                        values.append(len(stmt.evidence))
                    elif k == 'belief':
                        values.append(stmt.belief)
                    elif k == 'ag_count':
                        values.append(len(stmt.agent_list()))
                    else:
                        assert False, f"Value of k, {k}, should be impossible."
                stmt_metrics[sh] = array(values)
        else:
            # Iterate over the metric dict, and convert each entry to an array.
            for sh, metrics in metric_dict.items():
                stmt_metrics[sh] = array([metrics[k] for k in keys])

        new_cls = cls(keys, stmt_metrics, original_types)
        new_cls._set_metriks_to_stmt_metrics()
        new_cls.freeze()
        return new_cls

    def make_derivative_metriker(self):
        return self.__class__(self.__keys, self.__stmt_metrics,
                              self.__original_types)


class Metrik:
    def __init__(self, keys, stmt_metrics, original_types):
        self.keys = keys
        self.values = zeros(len(keys))
        self.stmt_metrics = stmt_metrics
        self.original_types = original_types
        self.__frozen = False

    @classmethod
    def from_array(cls, keys, arr, original_types, stmt_metrics=None):
        new_cls = cls(keys, stmt_metrics, original_types)
        new_cls.values = arr
        return new_cls

    def freeze(self):
        self.__frozen = True

    def include(self, stmt):
        if self.__frozen:
            raise RuntimeError("No longer adding more stmt data to Metrik.")
        if not isinstance(stmt, Statement):
            raise ValueError(f"Invalid type for addition to Metrik: "
                             f"{type(stmt)}. Must be a Statement.")

        h = stmt.get_hash()
        assert self.stmt_metrics and h in self.stmt_metrics
        self.values += self.stmt_metrics[h]

    def get_dict(self):
        return {key: value.astype(original_type)
                for key, value, original_type
                in zip(self.keys, self.values, self.original_types)}


def group_and_sort_statements(stmt_list, sort_by='default', metric_dict=None,
                              skip_grouping=False):
    """Group statements by type and arguments, and sort by prevalence.

    Parameters
    ----------
    stmt_list : list[Statement]
        A list of INDRA statements.
    sort_by : str or function
        If str, it indicates which parameter in metric_dict to sort by, or
        either 'belief' or 'ev_count' which are calculated from the statement
        objects themselves. The default, 'default', is mostly a sort by ev_count
        but also favors statements with fewer agents. Alternatively, you may
        give a function that takes a dict as its single argument, where that
        dict is an value of the `metric_dict`. If no metric dict is given,
        the argument to the function will recieve a dict with values for
        `belief` and `ev_count`.
    metric_dict : dict{int: dict}
        A dictionary of metrics that apply to each statement, keyed by hash.
        The value of `sort_by` should be a key into that dict, or should use the
        values of that dict (themselves a dict) as an argument if it is a
        function.
    skip_grouping : bool
        Optionally select to not group statements, and just sort them. Default
        is False.

    Returns
    -------
    sorted_groups : list[tuple]
        A list of tuples containing a sort key, the statement type, and a list
        of statements, also sorted by evidence count, for that key and type.
        The sort key contains a count of statements with those argument, the
        arguments (normalized strings), the count of statements with those
        arguments and type, and then the statement type.
    """
    relation_stmts = defaultdict(list)

    # Init the metrics.
    stmt_metrics = Metriker.from_stmt_list(stmt_list, metric_dict=metric_dict)
    relation_metrics = stmt_metrics.make_derivative_metriker()
    agent_pair_metrics = stmt_metrics.make_derivative_metriker()

    # Add up the grouped statements from the metrics.
    for rel_key, stmt in _get_relation_keyed_stmts(stmt_list):
        # Update the counts, and add key if needed.
        relation_stmts[rel_key].append(stmt)

        # Keep track of the evidence counts, source counts, and max belief for
        # this statement and the arguments.
        relation_metrics[rel_key].include(stmt)

        # Add up the counts for the arguments, pairwise for Complexes and
        # Conversions. This allows, for example, a complex between MEK, ERK,
        # and something else to lend weight to the interactions between MEK
        # and ERK.
        if rel_key[0] == 'Conversion':
            subj = rel_key[1]
            for obj in rel_key[2] + rel_key[3]:
                ag_pair_key = (subj, obj)
                agent_pair_metrics[ag_pair_key].include(stmt)
        else:
            ag_pair_key = rel_key[1:]
            agent_pair_metrics[ag_pair_key].include(stmt)

    # Stop filling these metrikers. No more "new" keys.
    relation_metrics.stop_filling()
    agent_pair_metrics.stop_filling()

    # Define the sort function.
    if isinstance(sort_by, str):
        def _sort_func(metric):
            assert isinstance(sort_by, str)
            if sort_by == 'default':
                return metric['ev_count'] + 1/(1 + metric['ag_count'])
            return metric[sort_by]
    else:
        _sort_func = sort_by

    # Sort the rows by count and agent names.
    def processed_rows(stmt_rows):
        for key, stmts in stmt_rows.items():
            verb = key[0]
            inps = key[1:]
            rel_m = relation_metrics[key]
            agp_m = agent_pair_metrics[inps]
            if verb == 'Complex' and rel_m['ev_count'] == agp_m['ev_count'] \
                    and len(inps) <= 2:
                if all([len(set(ag.name for ag in s.agent_list())) > 2
                        for s in stmts]):
                    continue
            agent_pair_sort_key = (_sort_func(agp_m.get_dict()), inps,
                                   _sort_func(rel_m.get_dict()), verb)

            def stmt_sorter(s):
                h = s.get_hash()
                metrics = stmt_metrics[h].get_dict()
                return _sort_func(metrics)

            stmts = sorted(stmts, key=stmt_sorter, reverse=True)
            yield agent_pair_sort_key, verb, stmts, \
                agent_pair_metrics[inps].get_dict(), \
                relation_metrics[key].get_dict()

    sorted_groups = sorted(processed_rows(relation_stmts),
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
    for key, s in _get_relation_keyed_stmts(stmts):
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
