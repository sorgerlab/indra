from collections import OrderedDict as _o
from datetime import timedelta


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


class QuantitativeState(Delta):
    """An object representing numerical value of something.

    Parameters
    ----------
    entity : str
        An entity to capture the quantity of.
    value : float or int
        Quantity of a unit (or range?)
    unit : str
        Measurement unit of value (e.g. absolute, daily, percentage, etc.)
    modifier : str
        Modifier to value (e.g. more than, at least, approximately, etc.)
    text : str
        Natural language text describing quantitative state.
    """
    def __init__(self, entity=None, value=None, unit=None, modifier=None,
                 text=None):
        self.entity = entity
        self.value = value
        self.unit = unit
        self.modifier = modifier
        self.text = text

    def equals(self, other):
        return (self.entity == other.entity and self.value == other.value and
                self.unit == other.unit and self.modifier == other.modifier and
                self.text == other.text)

    def to_json(self):
        json_dict = {'entity': self.entity if self.entity else None,
                     'value': self.value if self.value else None,
                     'unit': self.unit if self.unit else None,
                     'modifier': self.modifier if self.modifier else None,
                     'text': self.text if self.text else None}
        return json_dict

    @classmethod
    def from_json(cls, json_dict):
        entity = json_dict.get('entity')
        value = json_dict.get('value')
        unit = json_dict.get('unit')
        modifier = json_dict.get('modifier')
        text = json_dict.get('text')
        return cls(entity=entity, value=value, unit=unit, modifier=modifier,
                   text=text)

    def __str__(self):
        return "QuantitativeState(unit=%s, value=%d, text=%s)" % (self.unit,
                                                                  self.value,
                                                                  self.text)

    # Arithmetic operations
    def __add__(self, other, target_days=1):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for addition")
        values = self._standardize_units(other, target_unit='per_second')
        total_per_second = values[0] + values[1]
        total = self.from_seconds(total_per_second, target_days)
        return total

    def __sub__(self, other, target_days=1):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for subtraction")
        values = self._standardize_units(other, target_unit='per_second')
        diff_per_second = values[0] - values[1]
        diff = self.from_seconds(diff_per_second, target_days)
        return diff

    def __gt__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] > values[1])

    def __ge__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] >= values[1])

    def __lt__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] < values[1])

    def __le__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] <= values[1])

    def __eq__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] == values[1])

    def __ne__(self, other):
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        values = self._standardize_units(other, target_unit='per_second')
        return (values[0] != values[1])

    # Unit conversions
    def convert_unit(self, days, target_unit='per_second'):
        # convert between different units
        if target_unit == 'per_second':
            return self.value_per_second(self.value, days)

    def value_per_second(self, value, days):
        period = timedelta(days=days)
        seconds = period.total_seconds()
        per_second = value / seconds
        return per_second

    def from_seconds(self, value_per_second, days):
        period = timedelta(days=days)
        seconds = period.total_seconds()
        total_value = value_per_second * seconds
        return total_value

    def _standardize_units(self, other, target_unit='per_second'):
        values = []
        if self.unit != target_unit:
            days = self._get_days()
            self_new_value = self.convert_unit(days, target_unit)
            values.append(self_new_value)
        else:
            values.append(self.value)
        if other.unit != target_unit:
            days = other._get_days()
            other_new_value = other.convert_unit(days, target_unit)
            values.append(other_new_value)
        else:
            values.append(other.value)

    def _get_days(self):
        if self.unit == 'daily':
            return 1
        if self.unit == 'weekly':
            return 7
        if self.unit == 'monthly':
            return 30
        if self.unit == 'yearly':
            return 365
