__all__ = ['QueryResult', 'StatementQueryResult', 'AgentQueryResult']

from typing import Iterable as TypeIterable, Iterable, Union as TypeUnion

from indra.statements import stmts_from_json


class QueryResult(object):
    """The generic result of a query.

    This class standardizes the results of queries to the readonly database.

    Parameters
    ----------
    results : Iterable
        The results of the query keyed by unique IDs (mk_hash for PA Statements,
        IDs for Raw Statements, etc.)
    limit : int
        The limit that was applied to this query.
    offset_comp : int
        The next offset that would be appropriate if this is a paging query.
    evidence_counts : dict
        The count of evidence for each element.
    belief_scores : dict
        The belief score of each element.
    source_counts : dict
        The source counts for each element.
    query_json : dict
        A description of the query that was used.
    result_type : str
        The type of the result, e.g. 'agent', 'relation', etc.

    Attributes
    ----------
    results : Iterable
        The results of the query keyed by unique IDs (mk_hash for PA Statements,
        IDs for Raw Statements, etc.)
    limit : int
        The limit that was applied to this query.
    next_offset : int
        The next offset that would be appropriate if this is a paging query.
    evidence_counts : dict
        The count of evidence for each element.
    belief_scores : dict
        The belief score of each element.
    source_counts : dict
        The source counts for each element.
    query_json : dict
        A description of the query that was used.
    """
    def __init__(self, results: TypeIterable, limit: int, offset: int,
                 offset_comp: int, evidence_counts: dict, belief_scores: dict,
                 source_counts: dict, query_json: dict, result_type: str):
        if not isinstance(results, Iterable) or isinstance(results, str):
            raise ValueError("Input `results` is expected to be an iterable, "
                             "and not a string.")
        self.results = results
        self.evidence_counts = evidence_counts
        self.belief_scores = belief_scores
        self.total_evidence = sum(self.evidence_counts.values())
        if self.belief_scores:
            self.max_belief = max(self.belief_scores.values())
        else:
            self.max_belief = None
        self.source_counts = source_counts
        self.limit = limit
        self.offset = offset
        self.result_type = result_type
        self.offset_comp = offset_comp
        if limit is None or offset_comp < limit:
            self.next_offset = None
        else:
            self.next_offset = (0 if offset is None else offset) + offset_comp
        self.query_json = query_json

    @classmethod
    def empty(cls, empty_res, limit, offset, query_json, result_type):
        return cls(empty_res, limit, offset, 0, {}, {}, {}, query_json,
                   result_type)

    @classmethod
    def from_json(cls, json_dict) \
            -> TypeUnion['QueryResult', 'StatementQueryResult',
                         'AgentQueryResult']:
        # Build a StatementQueryResult or AgentQueryResult if appropriate
        if json_dict['result_type'] == 'statements':
            return StatementQueryResult.from_json(json_dict)
        elif json_dict['result_type'] == 'agents':
            return AgentQueryResult.from_json(json_dict)
        return cls._parse_json(json_dict)

    @classmethod
    def _parse_json(cls, json_dict):
        # Filter out some calculated values.
        next_offset = json_dict.pop('next_offset', None)
        total_evidence = json_dict.pop('total_evidence', None)

        # Build the class
        nc = cls(**json_dict)

        # Convert result keys into integers, if appropriate
        if nc.result_type in ['statements', 'interactions']:
            nc.results = {int(k): v for k, v in nc.results.items()}

        if nc.result_type in ['statements', 'interactions', 'hashes']:
            nc.evidence_counts = {int(k): v
                                  for k, v in nc.evidence_counts.items()}
            nc.belief_scores = {int(k): v for k, v in nc.belief_scores.items()}
            nc.source_counts = {int(k): v for k, v in nc.source_counts.items()}

        # Check calculated values.
        if nc.next_offset is None:
            nc.next_offset = next_offset
        else:
            assert nc.next_offset == next_offset, "Offsets don't match."
        assert nc.total_evidence == total_evidence,\
            "Evidence counts don't match."
        return nc

    def json(self) -> dict:
        """Return the JSON representation of the results."""
        if not isinstance(self.results, dict) \
                and not isinstance(self.results, list):
            json_results = list(self.results)
        elif isinstance(self.results, dict):
            json_results = {str(k): v for k, v in self.results.items()}
        else:
            json_results = self.results
        return {'results': json_results, 'limit': self.limit,
                'offset': self.offset, 'next_offset': self.next_offset,
                'query_json': self.query_json,
                'evidence_counts': self.evidence_counts,
                'belief_scores': self.belief_scores,
                'source_counts': self.source_counts,
                'total_evidence': self.total_evidence,
                'result_type': self.result_type,
                'offset_comp': self.offset_comp}


class StatementQueryResult(QueryResult):
    """The result of a query to retrieve Statements.

    This class encapsulates the results of a search for statements in the
    database. This standardizes the results of such searches.

    Parameters
    ----------
    results : dict
        The results of the query.
    limit : int
        The content limit that was used in making the query.
    offset : int
        The offset that was used in making the query.
    evidence_counts : dict
        The count of evidence available for each element.
    belief_scores : dict
        The score for each element.
    returned_evidence : int
        The count of evidence that was returned in this query.
    source_counts : dict
        The counts of evidence from each source for each element.
    query_json : dict
        The JSON representation of the query that was used.

    Attributes
    ----------
    results : dict
        The results of the query keyed by unique IDs (mk_hash for PA Statements,
        IDs for Raw Statements, etc.)
    limit : int
        The limit that was applied to this query.
    query_json : dict
        A description of the query that was used.
    """
    def __init__(self, results: dict, limit: int, offset: int,
                 evidence_counts: dict, belief_scores: dict,
                 returned_evidence: int, source_counts: dict, query_json: dict):
        super(StatementQueryResult, self).__init__(results, limit,
                                                   offset, len(results),
                                                   evidence_counts,
                                                   belief_scores, source_counts,
                                                   query_json, 'statements')
        self.returned_evidence = returned_evidence

    @classmethod
    def empty(cls, limit: int, offset: int, query_json: dict):
        return cls({}, limit, offset, {}, {}, 0, {}, query_json)

    def json(self) -> dict:
        """Get the JSON dump of the results."""
        json_dict = super(StatementQueryResult, self).json()
        json_dict.update({'returned_evidence': self.returned_evidence})
        return json_dict

    @classmethod
    def from_json(cls, json_dict):
        json_dict = json_dict.copy()
        result_type = json_dict.pop('result_type')
        json_dict.pop('offset_comp', None)
        if result_type != 'statements':
            raise ValueError(f'Invalid result type {result_type} for this '
                             f'result class {cls}')
        nc = super(StatementQueryResult, cls)._parse_json(json_dict)
        return nc

    def statements(self) -> list:
        """Get a list of Statements from the results."""
        assert isinstance(self.results, dict), "Results must be a dictionary."
        return stmts_from_json(list(self.results.values()))


class AgentQueryResult(QueryResult):
    """The result of a query for agent JSONs."""
    def __init__(self, results: dict, limit: int, offset: int, num_rows: int,
                 complexes_covered: set, evidence_counts: dict,
                 belief_scores: dict, source_counts: dict, query_json: dict):
        super(AgentQueryResult, self).__init__(results, limit, offset, num_rows,
                                               evidence_counts, belief_scores,
                                               source_counts, query_json,
                                               'agents')
        self.complexes_covered = complexes_covered

    @classmethod
    def empty(cls, limit, offset, query_json):
        return cls({}, limit, offset, 0, set(), {}, {}, {}, query_json)

    def json(self) -> dict:
        json_dict = super(AgentQueryResult, self).json()
        json_dict['complexes_covered'] = [str(h) for h in self.complexes_covered]
        return json_dict

    @classmethod
    def from_json(cls, json_dict):
        json_dict = json_dict.copy()
        result_type = json_dict.pop('result_type')
        if result_type != 'agents':
            raise ValueError(f'Invalid result type {result_type} for this '
                             f'result class {cls}')
        nc = super(AgentQueryResult, cls)._parse_json(json_dict)
        nc.complexes_covered = {int(h) for h in nc.complexes_covered}
