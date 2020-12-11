import logging
from collections import defaultdict
from itertools import permutations
from numpy import array, zeros, maximum, concatenate, append

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
        ags = s.agent_list()
        rel_key = None
        if verb == 'Complex':
            ag_ns = {name(ag) for ag in ags}
            if 1 < len(ag_ns) < 6:
                for pair in permutations(ag_ns, 2):
                    yield (verb,) + tuple(pair), tuple(pair), s
            if len(ag_ns) == 2:
                continue
            ag_key = tuple(sorted(ag_ns))
        elif verb == 'Conversion':
            subj = name(s.subj)
            objs_from = tuple(sorted({name(ag) for ag in s.obj_from}))
            objs_to = tuple(sorted({name(ag) for ag in s.obj_to}))
            for obj in objs_from:
                yield (verb, subj, objs_from, objs_to), (subj, obj), s
            for obj in objs_to:
                yield (verb, subj, objs_from, objs_to), (subj, obj), s
            ag_key = (subj, objs_from, objs_to)
        elif verb in ['ActiveForm', 'HasActivity']:
            ag_name = name(ags[0])
            ag_key = (ag_name,)
            rel_key = (verb, ag_name, s.activity,
                       s.is_active if verb == 'ActiveForm' else s.has_activity)
        elif verb == 'Influence':
            sns, sid = s.subj.concept.get_grounding()
            ons, oid = s.obj.concept.get_grounding()
            skey = s.subj.concept.name if not sid \
                else sid.split('/')[-1].replace('_', ' ')
            okey = s.obj.concept.name if not oid \
                else oid.split('/')[-1].replace('_', ' ')
            ag_key = (skey, okey)
        else:
            ag_key = tuple([name(ag) for ag in ags])

        # Set the default relation key.
        if rel_key is None:
            rel_key = (verb,) + ag_key

        # Yield the next (default) element.
        yield rel_key, ag_key, s


class StmtStat:
    """Abstraction of a metric applied to a set of statements.

    Parameters
    ----------
    name : str
        The label for this data (e.g. "ev_count" or "belief")
    data : dict{int: Number}
        The relevant statistics as a dict keyed by hash.
    data_type : type
        The type of the data (e.g. `int` or `float`).
    agg_class : type
        A subclass of BasicStats which defines how these statistics will be
        merged.
    """

    def __init__(self, name, data, data_type, agg_class):
        self.name = name
        self.data = data
        self.data_type = data_type
        self.agg_class = agg_class

    @classmethod
    def from_dicts(cls, dict_data, data_type, agg_class):
        """Generate a list of StmtStat's from a dict of dicts.

        Example Usage:
        >>> source_counts = {9623812756876: {'reach': 1, 'sparser': 2},
        >>>                  -39877587165298: {'reach': 3, 'sparser': 0}}
        >>> stmt_stats = StmtStat.from_dicts(source_counts, int, SumStats)

        Parameters
        ----------
        dict_data : dict{int: dict{str: Number}}
            A dictionary keyed by hash with dictionary elements, where each
            element gives a set of measurements for the statement labels as
            keys. A common example is `source_counts`.
        data_type : type
            The type of the data being given (e.g. `int` or `float`).
        agg_class : type
            A subclass of BasicStats which defines how these statistics will
            be merged (e.g. `SumStats`).
        """
        data_groups = defaultdict(dict)
        for h, data_dict in dict_data.items():
            for name, value in data_dict.items():
                data_groups[name][h] = value
            data_groups = dict(data_groups)

        classes = []
        for class_name, class_data in data_groups.items():
            classes.append(cls(class_name, class_data, data_type, agg_class))
        return classes

    @classmethod
    def from_stmts(cls, stmt_list, values=None):
        """Generate a list of StmtStat's from a list of stmts.

        The stats will include "ev_count", "belief", and "ag_count" by default,
        but a more limited selection may be specified using `values`.

        Example usage:
        >>> stmt_stats = StmtStat.from_stmts(stmt_list, ('ag_count', 'belief'))

        Parameters
        ----------
        stmt_list : list[Statement]
            A list of INDRA statements, from which basic stats will be derived.
        values : Optional[tuple(str)]
            A tuple of the names of the values to gather from the list of
            statements. For example, if you already have evidence counts, you
            might only want to gather belief and agent counts.
        """
        type_dict = {'ev_count': {'type': int, 'agg': SumStats},
                     'belief': {'type': float, 'agg': MaxStats},
                     'ag_count': {'type': int, 'agg': SumStats}}
        if values is None:
            values = tuple(type_dict.keys())

        # Iterate over statements, filling in values that may have been
        # missing.
        data = {k: {} for k in values}
        for stmt in stmt_list:
            sh = stmt.get_hash()
            if 'ev_count' in values:
                data['ev_count'][sh] = len(stmt.evidence)
            if 'belief' in values:
                data['belief'][sh] = stmt.belief
            if 'ag_count' in values:
                data['ag_count'][sh] = len(stmt.agent_list())

        # Create the objects.
        return [cls(k, d, type_dict[k]['type'], type_dict[k]['agg'])
                for k, d in data.items()]


class EvCount(StmtStat):
    def __init__(self, ev_counts):
        super(EvCount, self).__init__('ev_count', ev_counts, int, SumStats)


class Belief(StmtStat):
    def __init__(self, beliefs):
        super(Belief, self).__init__('belief', beliefs, float, MaxStats)


def source_count_list(source_counts):
    return StmtStat.from_dicts(source_counts, int, SumStats)


class StmtStatGather:
    """Gather metrics for items that are derived from statements.

    Example usage:
    >>> # Get ev_count, belief, and ag_count from a list of statements.
    >>> stmt_stats = StmtStat.from_stmts(stmt_list)
    >>>
    >>> # Add add another stat for a measure of relevance
    >>> stmt_stats.append(
    >>>     StmtStat('relevance', relevance_dict, float, AveStats)
    >>> )
    >>>
    >>> # Create the gatherer
    >>> StmtStatGather.from_stmt_stats(*stmt_stats)

    This class helps manage the accumulation of statistics for statements and
    statement-like objects, such as agent pairs. Working with StmtStatGroup and
    children of the BasicStats class, these tools make it easy to aggregate
    numerical measurements of statements with a great deal of flexibility.

    For example, you can sum up the evidence counts for statements that are part
    of an agent pair at the same time that you are remembering the maximum and
    average beliefs for that same corpus of statements. By defining your own
    child of BasicStats, specifically defining the operations that gather new
    data and finalize that data once all the statements are collected, you can
    utilize virtually any statistical methods for aggregating any metric for
    a Statement you might wish to use in sorting them.
    """

    def __init__(self, stat_groups):

        self.__stats = {}
        self.__started = False
        self.__finished = False

        # Check the groups and solidify them in more immutable types.
        hash_set = None
        self.__stmt_stats = {}
        rows = []
        for agg_class, info_dict in stat_groups.items():
            if hash_set is None:
                hash_set = set(info_dict['stats'].keys())
            else:
                if hash_set != set(info_dict['stats'].keys()):
                    raise ValueError(f"Stats from {info_dict['keys']} do "
                                     f"not cover the same corpora of hashes.")
            self.__stmt_stats[agg_class] = {
                'stats': {h: array(l) for h, l in info_dict['stats'].items()},
                'keys': tuple(info_dict['keys']),
                'types': tuple(info_dict['types'])
            }
            rows.extend(info_dict['keys'])

        self.__rows = tuple(rows)

    def add_stats(self, *stmt_stats):
        """Add more stats to the object."""
        if self.__started or self.__finished:
            raise RuntimeError("Cannot add stats after accumulation has "
                               "started or after it has finished.")

        for stat in stmt_stats:
            if not isinstance(stat, StmtStat):
                raise ValueError("All arguments must be StmtStat objects.")
            if stat.agg_class in self.__stmt_stats:
                self.__stmt_stats[stat.agg_class]['keys'] += (stat.name,)
                self.__stmt_stats[stat.agg_class]['types'] += (stat.data_type,)
                for h, v in stat.data.items():
                    old_arr = self.__stmt_stats[stat.agg_class]['stats'][h]
                    self.__stmt_stats[stat.agg_class]['stats'][h] = \
                        append(old_arr, v)
            else:
                self.__stmt_stats[stat.agg_class] = {
                    'stats': {h: array([v]) for h, v in stat.data.items()},
                    'keys': (stat.name,),
                    'types': (stat.data_type,)
                }
            self.__rows += (stat.name,)
        return

    def row_set(self):
        return set(self.__rows)

    @classmethod
    def from_stmt_stats(cls, *stmt_stats):
        """Create a stat gatherer from StmtStat objects."""

        # Organize the data into groups by aggregation class.
        stat_groups = defaultdict(lambda: {'stats': defaultdict(list),
                                           'keys': [], 'types': []})
        for stat in stmt_stats:
            if not isinstance(stat, StmtStat):
                raise ValueError("All arguments must be `StmtStat` object.")

            stat_groups[stat.agg_class]['keys'].append(stat.name)
            stat_groups[stat.agg_class]['types'].append(stat.data_type)
            for h, v in stat.data.items():
                stat_groups[stat.agg_class]['stats'][h].append(v)
        return cls(stat_groups)

    @classmethod
    def from_dicts(cls, ev_counts=None, beliefs=None, source_counts=None):
        """Init a stat gatherer from dicts keyed by hash."""
        stats = []
        if ev_counts:
            stats.append(EvCount(ev_counts))
        if beliefs:
            stats.append(Belief(beliefs))
        if source_counts:
            stats.extend(source_count_list(source_counts))
        return cls.from_stmt_stats(*stats)

    def __getitem__(self, key):
        if key not in self.__stats:
            if not self.__started:
                raise KeyError(f"Could not add key {key} before "
                               "accumulation started.")
            if not self.__finished:
                # Remember, this is passing REFERENCES to the stats dict.
                self.__stats[key] = StmtStatGroup(
                    agg_class(d['keys'], d['stats'], d['types'])
                    for agg_class, d in self.__stmt_stats.items()
                )
            else:
                raise KeyError(f"Key \"{key}\" not found! "
                               f"{self.__class__.__name__} is finished.")
        return self.__stats[key]

    def start(self):
        self.__started = True

    def finish(self):
        """Finish adding entries, new keys will be rejected."""
        self.__finished = True
        for stat_grp in self.__stats.values():
            stat_grp.finish()
        return

    def get_new_instance(self):
        """Create an instance to gather another level of data."""
        return self.__class__(self.__stmt_stats)

    def fill_from_stmt_stats(self):
        """Use the statements stats as stats and hashes as keys.

        This is used if you decide you just want to represent statements.
        """
        if self.__started or self.__finished:
            raise RuntimeError("Cannot fill from stats if accumulation has"
                               "already started or after it has finished.")

        # Gather stat rows from the stmt_stats.
        stat_rows = defaultdict(lambda: {'keys': tuple(), 'arr': array([]),
                                         'types': tuple()})
        for info_dict in self.__stmt_stats.values():
            for h, arr in info_dict['stats'].items():
                stat_rows[h]['keys'] += info_dict['keys']
                stat_rows[h]['arr'] = concatenate([stat_rows[h]['arr'], arr])
                stat_rows[h]['types'] += info_dict['types']
            stat_rows = dict(stat_rows)

        # Fill up the stats.
        for h, data in stat_rows.items():
            self.__stats[h] = BasicStats.from_array(data['keys'], data['arr'],
                                                    data['types'])

        # Mark as finished.
        self.finish()
        return


class StmtStatRow:
    """Define the API for elements of the stat"""

    def include(self, stmt):
        raise NotImplementedError()

    def get_dict(self):
        raise NotImplementedError()

    def finish(self):
        raise NotImplementedError()


class StmtStatGroup(StmtStatRow):
    """Implement the StmtStatRow API for a group of BasicStats children."""
    def __init__(self, stats):
        self.__stats = tuple(stats)
        self.__keymap = {k: stat for stat in self.__stats for k in stat.keys()}
        return

    def include(self, stmt):
        for stat in self.__stats:
            stat.include(stmt)

    def get_dict(self):
        return {k: v for stat in self.__stats
                for k, v in stat.get_dict().items()}

    def finish(self):
        for stat in self.__stats:
            stat.finish()

    def __getitem__(self, key):
        return self.__keymap[key][key]


class BasicStats(StmtStatRow):
    """Gathers measurements for a statement or similar entity.

    Parameters
    ----------
    keys : list[str]
        A dict keyed by aggregation method of lists of the names for the
        elements of data.
    stmt_metrics : dict{int: np.ndarray}
        A dictionary keyed by hash with each element a dict of arrays keyed
        by aggregation type.
    original_types : tuple(type)
        The type classes of each numerical value stored in the stmt_metrics
        dict, e.g. `(int, float, int)`.
    """

    def __init__(self, keys, stmt_metrics, original_types):
        self._keys = keys
        self._stmt_metrics = stmt_metrics
        self._original_types = original_types
        self._values = zeros(len(keys))
        self._count = 0
        self.__frozen = False

    @classmethod
    def from_array(cls, keys, arr, original_types, stmt_metrics=None):
        new_cls = cls(keys, stmt_metrics, original_types)
        new_cls._values = arr
        return new_cls

    def _finalize(self):
        return

    def finish(self):
        self._finalize()
        self.__frozen = True

    def include(self, stmt):
        if self.__frozen:
            raise RuntimeError("No longer adding more stmt data to BasicStats.")
        if not isinstance(stmt, Statement):
            raise ValueError(f"Invalid type for addition to BasicStats: "
                             f"{type(stmt)}. Must be a Statement.")

        h = stmt.get_hash()
        assert self._stmt_metrics and h in self._stmt_metrics
        self._merge(self._stmt_metrics[h])

    def _merge(self, metric_array):
        raise NotImplemented

    def __getitem__(self, item):
        if item not in self._keys:
            raise KeyError(f"Key '{item}' not found!")
        idx = self._keys.index(item)
        return self._values[idx].astype(self._original_types[idx])

    def keys(self):
        return self._keys[:]

    def get_dict(self):
        return {key: value.astype(original_type)
                for key, value, original_type
                in zip(self._keys, self._values, self._original_types)}


class SumStats(BasicStats):
    """A stats aggregator that executes a sum."""
    def _merge(self, metric_array):
        self._values += metric_array


class AveStats(BasicStats):
    """A stats aggregator averages the included statement metrics."""
    def _merge(self, metric_array):
        self._values += metric_array

    def _finalize(self):
        self._values = self._values / self._count


class MaxStats(BasicStats):
    """A stats aggregator that takes the max of statement metrics."""
    def _merge(self, metric_array):
        self._values = maximum(self._values, metric_array)


def _get_ag_name_set_len(stmt):
    return len(set(a.name if a else 'None' for a in stmt.agent_list()))


def group_and_sort_statements(stmt_list, sort_by='default', stmt_data=None,
                              grouping_level='agent-pair'):
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
    stmt_data : StmtStatGather
        A statement statistics gatherer loaded with data from the corpus of
        statements. If None, a new one will be formed with basic statics
        derived from the list of Statements itself.
    grouping_level : str
        The options are 'agent-pair', 'relation', and 'statement'. These
        correspond to grouping by agent pairs, agent and type relationships, and
        a flat list of statements. The default is 'agent-pair'.

    Returns
    -------
    sorted_groups : list[tuple]
        A list of tuples containing a sort key, the statement type, and a list
        of statements, also sorted by evidence count, for that key and type.
        The sort key contains a count of statements with those argument, the
        arguments (normalized strings), the count of statements with those
        arguments and type, and then the statement type.
    """
    # Validate the grouping level paramater.
    if grouping_level not in ['agent-pair', 'relation', 'statement']:
        raise ValueError(f"Invalid grouping level: \"{grouping_level}\".")

    # Init the stmt_data data gatherer.
    if stmt_data is None:
        stmt_stats = StmtStat.from_stmts(stmt_list)
        stmt_data = StmtStatGather.from_stmt_stats(*stmt_stats)
    missing_rows = {'ev_count', 'ag_count', 'belief'} - stmt_data.row_set()
    if missing_rows:
        stmt_stats = StmtStat.from_stmts(stmt_list, missing_rows)
        stmt_data.add_stats(*stmt_stats)
    stmt_data.fill_from_stmt_stats()

    # Define the sort function.
    if isinstance(sort_by, str):
        def _sort_func(metric):
            assert isinstance(sort_by, str)
            if sort_by == 'default':
                return metric['ev_count'] + 1/(1 + metric['ag_count'])
            return metric[sort_by]
    else:
        # Check that the sort function is a valid function.
        sample_dict = dict.fromkeys(stmt_data.row_set(), 0)
        try:
            n = sort_by(sample_dict)
            n < n
        except Exception as e:
            raise ValueError(f"Invalid sort function: {e}")

        # Assign the function.
        _sort_func = sort_by

    # Return the sorted statements, if that's all you want.
    if grouping_level == 'statement':
        def stmt_rows(stmts):
            for s in stmts:
                h = s.get_hash()
                metrics = stmt_data[h].get_dict()
                yield _sort_func(metrics), h, s, metrics
        return sorted(stmt_rows(stmt_list), key=lambda t: t[0], reverse=True)

    # Create gathering metrics from the statement data.
    relation_metrics = stmt_data.get_new_instance()
    relation_metrics.start()
    if grouping_level == 'agent-pair':
        agent_pair_metrics = stmt_data.get_new_instance()
        agent_pair_metrics.start()

    # Add up the grouped statements from the metrics.
    relation_stmts = defaultdict(list)
    for rel_key, ag_key, stmt in _get_relation_keyed_stmts(stmt_list):
        relation_metrics[rel_key].include(stmt)
        if grouping_level == 'agent-pair':
            relation_stmts[(ag_key, rel_key)].append(stmt)
            agent_pair_metrics[ag_key].include(stmt)
        else:
            relation_stmts[rel_key].append(stmt)

    # Stop filling these stat gatherers. No more "new" keys.
    relation_metrics.finish()
    if grouping_level == 'agent-pair':
        agent_pair_metrics.finish()

    # Sort the rows by count and agent names.
    def stmt_sorter(s):
        h = s.get_hash()
        metrics = stmt_data[h].get_dict()
        return _sort_func(metrics)

    if grouping_level == 'agent-pair':
        def processed_rows(stmt_rows):
            for (ag_key, rel_key), stmts in stmt_rows.items():
                verb = rel_key[0]
                rel_m = relation_metrics[rel_key]
                agp_m = agent_pair_metrics[ag_key]

                # Check if this statement is a type we ought to skip.
                is_abbrev_complex = verb == 'Complex' and len(ag_key) <= 2
                is_abbrev_conv = (verb == 'Conversion'
                                  and all(isinstance(e, str) for e in ag_key))
                is_abbrev_stmt = is_abbrev_conv or is_abbrev_complex
                all_ev_in_rel = rel_m['ev_count'] == agp_m['ev_count']
                all_have_many_ags = all(_get_ag_name_set_len(s) > 2
                                        for s in stmts)
                if is_abbrev_stmt and all_ev_in_rel and all_have_many_ags:
                    continue

                sort_key = (_sort_func(agp_m.get_dict()), str(ag_key),
                            _sort_func(rel_m.get_dict()), str(rel_key))

                stmts = sorted(stmts, key=stmt_sorter, reverse=True)

                yield sort_key, ag_key, rel_key, stmts, agp_m.get_dict(), \
                    rel_m.get_dict()
    else:  # If grouped by relation.
        def processed_rows(stmt_rows):
            for rel_key, stmts in stmt_rows.items():
                rel_m = relation_metrics[rel_key]
                sort_key = (_sort_func(rel_m.get_dict()), str(rel_key))

                yield sort_key, rel_key, stmts, rel_m.get_dict()

    return sorted(processed_rows(relation_stmts), key=lambda tpl: tpl[0],
                  reverse=True)


def make_stmt_from_relation_key(relation_key, agents=None):
    """Make a Statement from the relation key.

    Specifically, the sort key used by `group_and_sort_statements`.
    """

    def make_agent(name):
        if name == 'None' or name is None:
            return None
        return Agent(name)

    verb = relation_key[0]
    inps = relation_key[1:]
    StmtClass = get_statement_by_name(verb)
    if agents is None:
        agents = []
    if verb == 'Complex':
        agents.extend([make_agent(name) for name in inps])
        stmt = StmtClass(agents[:])
    elif verb == 'Conversion':
        if isinstance(inps[1], str):
            pass
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


def make_string_from_sort_key(rel_key):
    """Make a Statement string via EnglishAssembler from the sort key.

    Specifically, the sort key used by `group_and_sort_statements`.
    """
    stmt = make_stmt_from_relation_key(rel_key)
    return stmt_to_english(stmt)


def get_simplified_stmts(stmts):
    simple_stmts = []
    for rel_key, _, _ in _get_relation_keyed_stmts(stmts):
        simple_stmts.append(make_stmt_from_relation_key(rel_key))
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
