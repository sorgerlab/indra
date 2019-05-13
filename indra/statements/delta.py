class Delta(object):
    """The parent class of all delta types."""
    pass


class QualitativeDelta(Delta):
    """Qualitative delta defining an Event.

    Parameters
    ----------
    polarity : 1, -1 or None
        Polarity of an Event.
    adjectives : list[str]
        Adjectives describing an Event.
    """
    def __init__(self, polarity=None, adjectives=None):
        self.polarity = polarity
        self.adjectives = adjectives if adjectives else []

    def set_polarity(self, pol):
        self.polarity = pol

    def add_adjectives(self, adjectives):
        for adj in adjectives:
            self.adjectives.append(adj)

    def is_opposite(self, other):
        return ((self.polarity == 1 and other.polarity == -1) or
                (self.polarity == -1 and other.polarity == 1))
