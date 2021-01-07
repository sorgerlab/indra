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


def _get_relation_keyed_stmts(stmt_list, grouping_level='agent-pair'):
    def name(agent):
        return 'None' if agent is None else agent.name

    for s in stmt_list:
        # Create a key.
        verb = s.__class__.__name__
        ags = s.agent_list()
        rel_key = None
        if verb == 'Complex':
            ag_ns = {name(ag) for ag in ags}
            if grouping_level == 'agent-pair':
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
            if grouping_level == 'agent-pair':
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


def make_standard_stats(ev_counts=None, beliefs=None, source_counts=None):
    stats = []
    if ev_counts:
        stats.append(EvCount(ev_counts))
    if beliefs:
        stats.append(Belief(beliefs))
    if source_counts:
        stats.extend(source_count_list(source_counts))
    return stats


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
        new_stats = [s for s in stmt_stats if s.name not in self.row_set()]
        if not new_stats:
            return

        if self.__started or self.__finished:
            raise RuntimeError("Cannot add stats after accumulation has "
                               "started or after it has finished.")

        for stat in new_stats:
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
    def from_dicts(cls, **kwargs):
        """Init a stat gatherer from dicts keyed by hash."""
        stats = make_standard_stats(**kwargs)
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

    def is_finished(self):
        return self.__finished

    def is_started(self):
        return self.__started

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
        self.__stmt_hashes = set()
        self.__frozen = False
        self.__dict = None

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
        if h in self.__stmt_hashes:
            return
        assert self._stmt_metrics and h in self._stmt_metrics
        self._merge(self._stmt_metrics[h])
        self._count += 1
        self.__stmt_hashes.add(h)

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
        if not self.__frozen:
            raise RuntimeError("Cannot load source dict until frozen.")

        if self.__dict is None:
            self.__dict = {key: value.astype(original_type)
                           for key, value, original_type
                           in zip(self._keys, self._values,
                                  self._original_types)}
        return self.__dict


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


def group_and_sort_statements(stmt_list, sort_by='default', stmt_metrics=None,
                              grouping_level='agent-pair'):
    """Group statements by type and arguments, and sort by prevalence.

    Parameters
    ----------
    stmt_list : list[Statement]
        A list of INDRA statements.
    sort_by : str or function or None
        If str, it indicates which parameter to sort by, such as 'belief' or
        'ev_count', or 'ag_count'. Those are the default options because they
        can be derived from a list of statements, however if you give a custom
        `stmt_metrics`, you may use any of the parameters used to build it.
        The default, 'default', is mostly a sort by ev_count but also favors
        statements with fewer agents. Alternatively, you may give a function
        that takes a dict as its single argument, a dictionary of metrics. These
        metrics are determined by the contents of the `stmt_metrics` passed
        as an argument (see StmtStatGather for details), or else will contain
        the default metrics that can be derived from the statements themselves,
        namely `ev_count`, `belief`, and `ag_count`. The value may also
        be None, in which case the sort function will return the
        same value for all elements, and thus the original order of elements
        will be preserved. This could have strange effects when statements are
        grouped (i.e. when `grouping_level` is not 'statement'); such
        functionality is untested and we make no guarantee that it will work.
    stmt_metrics : StmtStatGather
        A statement statistics gatherer loaded with data from the corpus of
        statements. If None, a new one will be formed with basic statics
        derived from the list of Statements itself. See the docs for
        StmtStatGather for details on how to create one.
    grouping_level : str
        The options are 'agent-pair', 'relation', and 'statement'. These
        correspond to grouping by agent pairs, agent and type relationships, and
        a flat list of statements. The default is 'agent-pair'.

    Returns
    -------
    sorted_groups : list[tuple]
        A list of tuples of the form (sort_param, key, contents, metrics), where
        the sort param is whatever value was calculated to sort the results,
        the key is the unique key for the agent pair, relation, or statements,
        and the contents are either relations, statements, or statement JSON,
        depending on the level. This structure is recursive, so the each list
        of relations will also follow this structure, all the way down to
        the lowest level (statement JSON). The metrics a dict of the aggregated
        metrics for the entry (e.g. source counts, evidence counts, etc).
    """
    # Validate the grouping level parameter.
    if grouping_level not in ['agent-pair', 'relation', 'statement']:
        raise ValueError(f"Invalid grouping level: \"{grouping_level}\".")

    # Init the stmt_metrics data gatherer.
    if stmt_metrics is None:
        stmt_stats = StmtStat.from_stmts(stmt_list)
        stmt_metrics = StmtStatGather.from_stmt_stats(*stmt_stats)
    missing_rows = {'ev_count', 'ag_count', 'belief'} - stmt_metrics.row_set()
    if missing_rows:
        stmt_stats = StmtStat.from_stmts(stmt_list, missing_rows)
        stmt_metrics.add_stats(*stmt_stats)
    if not stmt_metrics.is_finished():
        stmt_metrics.fill_from_stmt_stats()

    # Define the sort function.
    if isinstance(sort_by, str):
        def _sort_func(metric):
            assert isinstance(sort_by, str)
            if sort_by == 'default':
                return metric['ev_count'] + 1/(1 + metric['ag_count'])
            return metric[sort_by]
    elif sort_by is None:
        def _sort_func(metric):
            return 0
    else:
        # Check that the sort function is a valid function.
        sample_dict = dict.fromkeys(stmt_metrics.row_set(), 0)
        try:
            n = sort_by(sample_dict)
            n < n
        except Exception as e:
            raise ValueError(f"Invalid sort function: {e}")

        # Assign the function.
        _sort_func = sort_by

    # Return the sorted statements, if that's all you want.
    def iter_rows(rows, *metric_dicts):
        assert metric_dicts
        for key, contents in rows:
            metrics = metric_dicts[0][key].get_dict()
            if len(metric_dicts) > 1:
                if isinstance(contents, dict):
                    contents = contents.items()
                contents = sorted_rows(contents, *metric_dicts[1:])
            yield (_sort_func(metrics), str(key)) if sort_by else 0, \
                  key, contents, metrics

    def sorted_rows(rows, *metric_dicts):
        return sorted(iter_rows(rows, *metric_dicts), key=lambda t: t[0],
                      reverse=True)

    if grouping_level == 'statement':
        stmt_rows = ((s.get_hash(), s) for s in stmt_list)
        return sorted_rows(stmt_rows, stmt_metrics)

    # Create gathering metrics from the statement data.
    relation_metrics = stmt_metrics.get_new_instance()
    relation_metrics.start()
    if grouping_level == 'agent-pair':
        agent_pair_metrics = stmt_metrics.get_new_instance()
        agent_pair_metrics.start()

    # Add up the grouped statements from the metrics.
    if grouping_level == 'relation':
        grouped_stmts = defaultdict(list)
    else:
        grouped_stmts = defaultdict(lambda: defaultdict(list))
    for rel_key, ag_key, stmt in _get_relation_keyed_stmts(stmt_list,
                                                           grouping_level):
        relation_metrics[rel_key].include(stmt)
        if grouping_level == 'agent-pair':
            grouped_stmts[ag_key][rel_key].append((stmt.get_hash(), stmt))
            agent_pair_metrics[ag_key].include(stmt)
        else:
            grouped_stmts[rel_key].append((stmt.get_hash(), stmt))

    # Stop filling these stat gatherers. No more "new" keys.
    relation_metrics.finish()
    if grouping_level == 'agent-pair':
        agent_pair_metrics.finish()

    # Sort the rows by count and agent names.
    if grouping_level == 'relation':
        return sorted_rows(grouped_stmts.items(), relation_metrics,
                           stmt_metrics)

    return sorted_rows(grouped_stmts.items(), agent_pair_metrics,
                       relation_metrics, stmt_metrics)


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


def make_string_from_relation_key(rel_key):
    """Make a Statement string via EnglishAssembler from the relation key.

    Specifically, the key used by `group_and_sort_statements` for contents
    grouped at the relation level.
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
