from collections import OrderedDict as _o


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

    def equals(self, other):
        return (self.polarity == other.polarity and
                set(self.adjectives) == set(other.adjectives))

    def to_json(self):
        json_dict = _o({'polarity': self.polarity})
        if self.adjectives:
            json_dict['adjectives'] = self.adjectives
        return json_dict

    @classmethod
    def from_json(cls, json_dict):
        polarity = json_dict.get('polarity')
        adjectives = json_dict.get('adjectives')
        delta = QualitativeDelta(polarity=polarity, adjectives=adjectives)
        return delta

    def __str__(self):
        return "%s(polarity=%d, adjectives=%s)" % (type(self).__name__,
                                                   self.polarity,
                                                   self.adjectives)
