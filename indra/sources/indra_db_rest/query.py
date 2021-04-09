"""
The Query architecture allows the construction of arbitrary queries for content
from the INDRA Database.

Specifically, the query constructed using this language of classes is converted
into sophisticated and optimized SQL. Different classes represent different
types of constraint and are named as much as possible to fit together when
spoken aloud in English. For example:

>>> HasAgent("MEK") & HasAgent("ERK") & HasType(["Phosphorylation"])

Will find any Statement that has an agent MEK and an agent ERK and has the type
phosphorylation. Broadly, querie classes can be broken into 3 types: queries on
the meaning of a Statement, queries on the provenance of a Statement, and
queries that combine groups of queries.

Meaning of a Statement:

- :py:class:`HasAgent`
- :py:class:`HasType`
- :py:class:`HasNumAgents`

Provenence of a Statement:

- :py:class:`HasReadings`
- :py:class:`HasDatabases`
- :py:class:`HasSources`
- :py:class:`HasOnlySource`
- :py:class:`FromPapers`
- :py:class:`FromMeshIds`
- :py:class:`HasNumEvidence`
- :py:class:`HasEvidenceBound`

Combine Queriers:

- :py:class:`And`
- :py:class:`Or`

There is also the special class, the :py:class:`EmptyQuery` which is useful when
programmatically building a query.

In practice you will likely not use :py:class:`And` or :py:class:`Or` very often
but will instead make use of the overloaded `&` and `|` operators. In addition
you can invert a query -- essentially ask for Statements that do _not_ meet
certain criteria, e.g. "not has readings". This can be accomplished with the
overloaded `~` operator, e.g. `~HasReadings()`.

The query class works by representing a and producing a particular JSON
structure which is recognized by the INDRA Database REST service, where it is
translated into a similar but more sophisticated Query language used by the
Readonly Database client. The Query class implements the basic methods used to
communicate with the REST Service in this way.

Examples
--------

First an example of the typical usage of a query object:

>>> from indra.sources.indra_db_rest.api import get_statements_from_query
>>> from indra.sources.indra_db_rest.query import *
>>>
>>> q = HasAgent('MEK') | HasAgent('MAP2K1') & HasDatabases()
>>> p = get_statements_from_query(q)
>>> p.statements
[Activation(MEK(), ERK()),
 Phosphorylation(MEK(), ERK()),
 Activation(MAP2K1(), ERK()),
 Activation(RAF1(), MEK()),
 Phosphorylation(RAF1(), MEK()),
 Phosphorylation(MAP2K1(), ERK()),
 Activation(BRAF(), MEK()),
 Inhibition(2-(2-amino-3-methoxyphenyl)chromen-4-one(), MEK()),
 Activation(MAP2K1(), MAPK1()),
 Activation(MAP2K1(), MAPK3()),
 Phosphorylation(MAP2K1(), MAPK1()),
 Phosphorylation(BRAF(), MEK()),
 Activation(MEK(), MAPK1()),
 Complex(BRAF(), MAP2K1()),
 Phosphorylation(MAP2K1(), MAPK3()),
 Activation(MEK(), MAPK3()),
 Complex(MAP2K1(), RAF1()),
 Activation(RAF1(), MAP2K1()),
 Inhibition(trametinib(), MEK()),
 Phosphorylation(MEK(), MAPK3()),
 Complex(MAP2K1(), MAPK1()),
 Phosphorylation(MEK(), MAPK1()),
 Inhibition(selumetinib(), MEK()),
 Phosphorylation(PAK1(), MAP2K1(), S, 298)]
>>>
>>> q = HasAgent('MEK') & HasAgent('ERK') & HasEvidenceBound(["> 10"])
>>> p = get_statements_from_query(q)
>>> p.statements
[Activation(MEK(), ERK()),
 Phosphorylation(MEK(), ERK()),
 Complex(ERK(), MEK()),
 Inhibition(MEK(), ERK()),
 Dephosphorylation(MEK(), ERK()),
 Complex(ERK(), MEK(), RAF()),
 Phosphorylation(MEK(), ERK(), T),
 Phosphorylation(MEK(), ERK(), Y),
 Activation(MEK(), ERK(mods: (phosphorylation))),
 IncreaseAmount(MEK(), ERK())]

For more details on the usage of :py:function:`get_statments_form_query` see the
:py:module:`indra.sources.indra_db_rest.api` documentation.

Lastly, a tour demonstrating the different utilities of a query object:

>>> from indra.sources.indra_db_rest.query import *
>>> q = HasAgent('MEK', namespace='FPLX') & ~HasAgent('ERK', namespace='FPLX')
>>>
>>> # This is the JSON sent to the server.
>>> q.to_simple_json()
{'class': 'And',
 'constraint': {'queries': [{'class': 'HasAgent',
    'constraint': {'agent_id': 'MEK',
     'namespace': 'FPLX',
     'role': None,
     'agent_num': None},
    'inverted': False},
   {'class': 'HasAgent',
    'constraint': {'agent_id': 'ERK',
     'namespace': 'FPLX',
     'role': None,
     'agent_num': None},
    'inverted': True}]},
 'inverted': False}
>>>
>>> # The more "true" representation of the JSON in this case looks very similar
>>> q.get_query_json()
{'class': 'Intersection',
 'constraint': {'query_list': [{'class': 'HasAgent',
    'constraint': {'_regularized_id': 'MEK',
     'agent_id': 'MEK',
     'agent_num': None,
     'namespace': 'FPLX',
     'role': None},
    'inverted': False},
   {'class': 'HasAgent',
    'constraint': {'_regularized_id': 'ERK',
     'agent_id': 'ERK',
     'agent_num': None,
     'namespace': 'FPLX',
     'role': None},
    'inverted': True}]},
 'inverted': False}
>>>
>>> print("I am finding statements that", q.get_query_english())
I am finding statements that do not have an agent where FPLX=ERK and have an
agent where FPLX=MEK
"""

__all__ = ['Query', 'And', 'Or', 'HasAgent', 'FromMeshIds', 'HasHash',
           'HasSources', 'HasOnlySource', 'HasReadings', 'HasDatabases',
           'HasType', 'HasNumAgents', 'HasNumEvidence', 'HasEvidenceBound',
           'FromPapers', 'EmptyQuery']

from typing import Iterable, Tuple, Union

from indra.sources.indra_db_rest.query_results import QueryResult
from indra.sources.indra_db_rest.util import make_db_rest_request, jsonify_args


class Query:
    """The parent of all query objects."""
    def __init__(self):
        self._inverted = False
        self.__compiled_json = None
        self.__compiled_str = None

    # Here are defined some other functions to get info from the server.

    def get(self, result_type, limit=None, sort_by=None, offset=None,
            timeout=None, n_tries=2, api_key=None, **other_params):
        """Get results from the API of the given type.

        Parameters
        ----------
        result_type : str
            The options are 'statements', 'interactions', 'relations', 'agents',
            and 'hashes', indicating the type of result you want.
        limit : Optional[int]
            The maximum number of statements you want to try and retrieve. The
            server will by default limit the results, and any value exceeding
            that limit will be "overruled".
        sort_by : Optional[str]
            The value can be 'default', 'ev_count', or 'belief'.
        offset : Optional[int]
            The offset of the query to begin at.
        timeout : Optional[int]
            The number of seconds to wait for the request to return before
            giving up. This timeout is applied to each try separately.
        n_tries : Optional[int]
            The number of times to retry the request before giving up. Each try
            will have `timeout` seconds to complete before it gives up.
        api_key : str or None
            Override or use in place of the API key given in the INDRA config
            file.

        Other Parameters
        ----------------
        'statements'
        filter_ev : bool
            Indicate whether evidence should have the same filters applied as
            the statements themselves, where appropriate (e.g. in the case of a
            filter by paper).
        ev_limit : int
            Limit the number of evidence returned per Statement.

        'relations' and 'agents'
        with_hashes : bool
            Choose whether the hashes for each Statement be included along with
            each grouped heading.

        'agents'
        complexes_covered : list[int]
            A list (or set) of complexes that have already come up in the agent
            groups returned. This prevents duplication.
        """
        simple = self.__compiled_json is None
        if simple:
            query_json = self.to_simple_json()
        else:
            query_json = self.__compiled_json
        resp = make_db_rest_request('post', f'query/{result_type}',
                                    data={'query': query_json,
                                          'kwargs': jsonify_args(other_params)},
                                    params=dict(limit=limit, sort_by=sort_by,
                                                offset=offset, simple=simple),
                                    timeout=timeout, tries=n_tries,
                                    api_key=api_key)
        resp_json = resp.json()
        self.__compiled_json = resp_json['query_json']
        self.__compiled_str = None
        return QueryResult.from_json(resp_json)

    def get_query_json(self):
        """Generate a compiled JSON rep of the query on the server."""
        if not self.__compiled_json:
            resp = make_db_rest_request('post', 'compile/json',
                                        data=self.to_simple_json())
            self.__compiled_json = resp.json()
            self.__compiled_str = None
        return self.__compiled_json

    def get_query_english(self):
        """Get the string representation of the query."""
        if self.__compiled_str is None:
            if self.__compiled_json is None:
                query_json = self.to_simple_json()
                simple = True
            else:
                query_json = self.__compiled_json
                simple = False
            resp = make_db_rest_request('post', 'compile/string',
                                        data=query_json,
                                        params=dict(simple=simple))
            self.__compiled_str = resp.content.decode('utf-8')
        return self.__compiled_str

    # Local (and largely internal) tools:

    def copy(self):
        """Make a copy of the query."""
        cp = self._copy()
        cp._inverted = self._inverted
        return cp

    def _copy(self):
        raise NotImplementedError()

    def to_simple_json(self) -> dict:
        """Generate the JSON from the object rep."""
        return {'class': self.__class__.__name__,
                'constraint': self.get_constraint_dict(),
                'inverted': self._inverted}

    def get_constraint_dict(self) -> dict:
        raise NotImplementedError()

    def invert(self):
        return self.__invert__()

    # Define the operator overloads.

    def __and__(self, other):
        if isinstance(other, EmptyQuery):
            return self.copy()
        return And([self.copy(), other.copy()])

    def __or__(self, other):
        if isinstance(other, EmptyQuery):
            return self.copy()
        return Or([self.copy(), other.copy()])

    def __invert__(self):
        inv = self.copy()
        inv._inverted = not self._inverted

        return inv

    def __repr__(self):
        inv = '~' if self._inverted else ''
        args = ', '.join(f'{key}="{value}"' if isinstance(value, str)
                         else f'{key}={value}'
                         for key, value in self.get_constraint_dict().items())
        return f"{inv}{self.__class__.__name__}({args})"


class And(Query):
    """The intersection of two queries.

    This are generally generated from the use of &, e.g.
    q_and = HashAgent('MEK') & HasAgent('ERK').
    """

    def __init__(self, queries: list):
        self.queries = queries
        super(And, self).__init__()

    def _copy(self):
        return And([q.copy() for q in self.queries])

    def get_constraint_dict(self) -> dict:
        return {'queries': [q.to_simple_json() for q in self.queries]}

    def __repr__(self):
        q_strings = [repr(q) for q in self.queries]
        s = ' & '.join(q_strings)
        if self._inverted:
            s = f'~({s})'
        return s

    def __and__(self, other):
        if isinstance(other, And):
            other_queries = other.queries
        else:
            other_queries = [other]
        return And([q.copy() for q in (self.queries + other_queries)])


class Or(Query):
    """The union of two queries.

    These are generally generate from the use of '|', e.g.
    q_or = HasOnlySource('reach') | HasOnlySource('medscan').
    """

    def __init__(self, queries: list):
        self.queries = queries
        super(Or, self).__init__()

    def _copy(self):
        return Or([q.copy() for q in self.queries])

    def get_constraint_dict(self) -> dict:
        return {'queries': [q.to_simple_json() for q in self.queries]}

    def __repr__(self):
        q_strings = [repr(q) for q in self.queries]
        s = ' | '.join(q_strings)
        if self._inverted:
            s = f'~({s})'
        return s

    def __or__(self, other):
        if isinstance(other, Or):
            other_queries = other.queries
        else:
            other_queries = [other]
        return Or([q.copy() for q in (self.queries + other_queries)])


class EmptyQuery(Query):
    """A query that is empty."""

    def _copy(self):
        return EmptyQuery()

    def __and__(self, other):
        return other

    def __or__(self, other):
        return other

    def get_constraint_dict(self) -> dict:
        return {}


class HasOnlySource(Query):
    """Find Statements that come exclusively from one source.

    For example, find statements that come only from sparser.

    Parameters
    ----------
    only_source : str
        The only source that spawned the statement, e.g. signor, or reach.
    """

    def __init__(self, only_source):
        self.only_source = only_source
        super(HasOnlySource, self).__init__()

    def _copy(self):
        return HasOnlySource(self.only_source)

    def get_constraint_dict(self) -> dict:
        return {'only_source': self.only_source}


class HasSources(Query):
    """Find Statements with support from the given list of sources.

    For example, find Statements that have support from both medscan and reach.

    Parameters
    ----------
    sources : list or set or tuple
        A collection of strings, each string the canonical name for a source.
        The result will include statements that have evidence from ALL sources
        that you include.
    """

    def __init__(self, sources):
        self.sources = tuple(set(sources))
        super(HasSources, self).__init__()

    def _copy(self):
        return HasSources(self.sources[:])

    def get_constraint_dict(self) -> dict:
        return {'sources': self.sources}


class HasReadings(Query):
    """Find Statements with support from readings."""

    def _copy(self):
        return HasReadings()

    def get_constraint_dict(self) -> dict:
        return {}


class HasDatabases(Query):
    """Find Statements with support from Databases."""

    def _copy(self):
        return HasDatabases()

    def get_constraint_dict(self) -> dict:
        return {}


class HasHash(Query):
    """Find Statements whose hash is contained in the given list.

    Parameters
    ----------
    stmt_hashes : list or set or tuple
        A collection of integers, where each integer is a shallow matches key
        hash of a Statement (frequently simply called "mk_hash" or "hash")
    """

    def __init__(self, stmt_hashes):
        self.stmt_hashes = stmt_hashes
        super(HasHash, self).__init__()

    def _copy(self):
        return HasHash(self.stmt_hashes)

    def get_constraint_dict(self) -> dict:
        return {'stmt_hashes': self.stmt_hashes}


class HasAgent(Query):
    """Find Statements with the given agent in the given position.

    Parameters
    ----------
    agent_id : str
        The ID string naming the agent, for example 'ERK' (FPLX or NAME) or
        'plx' (TEXT), and so on.
    namespace : str
        (optional) By default, this is AUTO, indicating GILDA will be used to
        to try and guess the proper namespace and agent ID. Other options
        include NAME (the canonical name of the agent), FPLX (FamPlex), CHEBI,
        CHEMBL, HGNC, UP (UniProt), TEXT (for raw text mentions), and many more.
    role : str or None
        (optional) None by default. Options are "SUBJECT", "OBJECT", or "OTHER".
    agent_num : int or None
        (optional) None by default. The regularized position of the agent in the
        Statement's list of agents.
    """

    def __init__(self, agent_id, namespace='NAME', role=None, agent_num=None):
        self.agent_id = agent_id
        self.namespace = namespace
        self.role = role
        self.agent_num = agent_num
        super(HasAgent, self).__init__()

    def _copy(self):
        return HasAgent(self.agent_id, self.namespace, self.role,
                        self.agent_num)

    def get_constraint_dict(self) -> dict:
        return {'agent_id': self.agent_id, 'namespace': self.namespace,
                'role': self.role, 'agent_num': self.agent_num}


class FromPapers(Query):
    """Get Statements that came from a given list of papers.

    Parameters
    ----------
    paper_list : list[(<id_type>, <paper_id>)]
        A list of tuples, where each tuple indicates and id-type (e.g. 'pmid')
        and an id value for a particular paper.
    """

    def __init__(self, paper_list):
        self.paper_list = paper_list
        super(FromPapers, self).__init__()

    def _copy(self):
        return FromPapers(self.paper_list)

    def get_constraint_dict(self) -> dict:
        return {'paper_list': self.paper_list}


class FromMeshIds(Query):
    """Get stmts that came from papers annotated with the given Mesh Ids.

    Parameters
    ----------
    mesh_ids : list
        A canonical MeSH ID, of the "C" or "D" variety, e.g. "D000135".
    """

    def __init__(self, mesh_ids):
        self.mesh_ids = mesh_ids
        super(FromMeshIds, self).__init__()

    def _copy(self):
        return FromMeshIds(self.mesh_ids)

    def get_constraint_dict(self) -> dict:
        return {'mesh_ids': self.mesh_ids}


class HasNumAgents(Query):
    """Get Statements with the given number of agents.

    For example, `HasNumAgents([1,3,4])` will return agents with either 2,
    3, or 4 agents (the latter two mostly being complexes).

    Parameters
    ----------
    agent_nums : tuple
        A list of integers, each indicating a number of agents.
    """

    def __init__(self, agent_nums):
        self.agent_nums = agent_nums
        super(HasNumAgents, self).__init__()

    def _copy(self):
        return HasNumAgents(self.agent_nums)

    def get_constraint_dict(self) -> dict:
        return {'agent_nums': self.agent_nums}


class HasNumEvidence(Query):
    """Get Statements with the given number of evidence.

    For example, HasNumEvidence([2,3,4]) will return Statements that have
    either 2, 3, or 4 evidence.

    Parameters
    ----------
    evidence_nums :
        A list of numbers greater than 0, each indicating a number of evidence.
    """

    def __init__(self, evidence_nums: Tuple[Union[int, str]]):
        self.evidence_nums = evidence_nums
        super(HasNumEvidence, self).__init__()

    def _copy(self):
        return HasNumEvidence(self.evidence_nums)

    def get_constraint_dict(self) -> dict:
        return {'evidence_nums': self.evidence_nums}


class HasEvidenceBound(Query):
    """Get Statements with given bounds on their evidence count.

    For example, HasEvidenceBound(["< 10", ">= 5"]) will return Statements with
    less than 10 and as many or more than 5 evidence.

    Parameters
    ----------
    evidence_bounds :
        An iterable (e.g. list) of strings such as "< 2" or ">= 4". The argument
        of the inequality must be a natural number (0, 1, 2, ...) and the
        inequality operation must be one of: <, >, <=, >=, ==, !=.
    """

    def __init__(self, evidence_bounds: Union[Iterable[str], str]):
        if isinstance(evidence_bounds, str):
            evidence_bounds = [evidence_bounds]
        self.evidence_bounds = list(evidence_bounds)
        super(HasEvidenceBound, self).__init__()

    def _copy(self):
        return HasEvidenceBound(self.evidence_bounds)

    def get_constraint_dict(self) -> dict:
        return {'evidence_bounds': self.evidence_bounds}


class HasType(Query):
    """Get Statements with the given type.


    For example, you can find Statements that are Phosphorylations or
    Activations, or you could find all subclasses of RegulateActivity.

    Parameters
    ----------
    stmt_types : set or list or tuple
        A collection of Strings, where each string is a class name for a type
        of Statement. Spelling and capitalization are necessary.
    include_subclasses : bool
        (optional) default is False. If True, each Statement type given in the
        list will be expanded to include all of its sub classes.
    """

    def __init__(self, stmt_types, include_subclasses=False):
        if isinstance(stmt_types, str):
            stmt_types = [stmt_types]
        self.stmt_types = stmt_types
        self.include_subclasses = include_subclasses
        super(HasType, self).__init__()

    def _copy(self):
        return HasType(self.stmt_types, self.include_subclasses)

    def get_constraint_dict(self) -> dict:
        return {'stmt_types': self.stmt_types,
                'include_subclasses': self.include_subclasses}
