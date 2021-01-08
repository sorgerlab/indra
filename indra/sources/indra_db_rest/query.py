class Query:
    """The parent of all query objects."""
    def __init__(self):
        self._inverted = False

    # Here are defined some other functions to get info from the server.

    def compile(self):
        """Generate a compiled JSON rep of the query on the server."""

    # Local (and largely internal) tools:

    def copy(self):
        """Make a copy of the query."""
        cp = self._copy()
        cp._inverted = self._inverted
        return cp

    def _copy(self):
        raise NotImplementedError()

    def to_json(self) -> dict:
        """Generate the JSON from the object rep."""
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


class And(Query):
    """The intersection of two queries."""

    def __init__(self, queries: list):
        self.queries = queries
        super(And, self).__init__()

    def _copy(self):
        return And([q.copy() for q in self.queries])

    def to_json(self) -> dict:
        return {'type': 'And',
                'constraints': {'queries': [q.to_json() for q in self.queries]},
                'inverted': self._inverted}


class Or(Query):
    """The union of two queries."""

    def __init__(self, queries: list):
        self.queries = queries
        super(Or, self).__init__()

    def _copy(self):
        return Or([q.copy() for q in self.queries])

    def to_json(self) -> dict:
        return {'type': 'Or',
                'constraints': {'queries': [q.to_json() for q in self.queries]},
                'inverted': self._inverted}


class EmptyQuery(Query):
    """A query that is empty."""

    def _copy(self):
        return EmptyQuery()

    def __and__(self, other):
        return other

    def __or__(self, other):
        return other

    def to_json(self) -> dict:
        return {}


class HasOnlySource(Query):
    """Find Statements that come exclusively from one source."""

    def __init__(self, only_source):
        self.only_source = only_source
        super(HasOnlySource, self).__init__()

    def _copy(self):
        return HasOnlySource(self.only_source)


class HasSources(Query):
    """Find Statements with support from the given list of sources."""


class HasReadings(Query):
    """Find Statements with support from readings."""


class HasDatabases(Query):
    """Find Statements with support from Databases."""


class HasHash(Query):
    """Find Statements whose hash is contained in the given list."""


class HasAgent(Query):
    """Find Statements with the given agent in the given position."""


class FromPapers(Query):
    """Get Statements that came from a given list of papers."""


class FromMeshIds(Query):
    """Get Statements that came from papers annotated with the given Mesh Ids."""


class HasNumAgents(Query):
    """Get Statements with the given number of agents."""


class HasNumEvidence(Query):
    """Get Statements with the given number of evidence."""


class HasType(Query):
    """Get Statements with the given type."""
