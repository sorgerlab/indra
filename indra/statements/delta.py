import logging
from collections import OrderedDict as _o
from datetime import timedelta


__all__ = ['Delta', 'QualitativeDelta', 'QuantitativeState']


logger = logging.getLogger(__name__)


class Delta(object):
    """The parent class of all delta types."""
    @classmethod
    def from_json(cls, json_dict):
        delta_type = json_dict.get('type')
        if delta_type == 'qualitative':
            return QualitativeDelta.from_json(json_dict)
        elif delta_type == 'quantitative':
            return QuantitativeState.from_json(json_dict)
        else:
            raise ValueError('Unknown delta type %s' % delta_type)

    def refinement_of(self, other):
        if type(self) != type(other):
            return False
        if self.polarity and not other.polarity:
            return True
        if self.polarity == other.polarity:
            return True

    def set_polarity(self, pol):
        self.polarity = pol


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
        json_dict = _o({'type': 'qualitative', 'polarity': self.polarity})
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
        return "%s(polarity=%s, adjectives=%s)" % (type(self).__name__,
                                                   str(self.polarity),
                                                   self.adjectives)

    def __repr__(self):
        return self.__str__()


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
    polarity : 1, -1 or None
        Polarity of an Event.
    """
    def __init__(self, entity=None, value=None, unit=None, modifier=None,
                 text=None, polarity=None):
        self.entity = entity
        self.value = value if value else None
        self.unit = unit
        self.modifier = modifier
        self.text = text
        self.polarity = polarity

    def equals(self, other):
        return (self.entity == other.entity and self.value == other.value and
                self.unit == other.unit and self.modifier == other.modifier and
                self.text == other.text and self.polarity == other.polarity)

    def refinement_of(self, other):
        ref = super().refinement_of(other)
        if not ref:
            return False
        if self.entity != other.entity:
            return False
        if not other.value:
            return True
        elif not self.value:
            return False
        if self.value == other.value and self.unit == other.unit:
            return True
        else:
            return False

    def to_json(self):
        json_dict = {'type': 'quantitative',
                     'entity': self.entity if self.entity else None,
                     'value': self.value if self.value else None,
                     'unit': self.unit if self.unit else None,
                     'modifier': self.modifier if self.modifier else None,
                     'text': self.text if self.text else None,
                     'polarity': self.polarity if self.polarity else None}
        return json_dict

    @classmethod
    def from_json(cls, json_dict):
        entity = json_dict.get('entity')
        value = json_dict.get('value')
        unit = json_dict.get('unit')
        modifier = json_dict.get('modifier')
        text = json_dict.get('text')
        polarity = json_dict.get('polarity')
        return cls(entity=entity, value=value, unit=unit, modifier=modifier,
                   text=text, polarity=polarity)

    def __str__(self):
        return ("QuantitativeState(entity=%s, value=%s, unit=%s, modifier=%s,"
                " text=%s, polarity=%s)" % (
                    self.entity, str(self.value), self.unit,
                    self.modifier, self.text, str(self.polarity)))

    def __repr__(self):
        return self.__str__()

    # Arithmetic operations
    def __add__(self, other):
        # Return value and unit
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for addition")
        # If original units are the same, return the result in same unit
        if self.unit == other.unit:
            return ((self.value + other.value), self.unit)
        # Otherwise, return result per second
        values = self._standardize_units(other, target_unit='second')
        total_per_second = values[0] + values[1]
        # result is per second
        return (total_per_second, 'second')

    def __sub__(self, other):
        # Return value and unit
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for subtraction")
        # If original units are the same, return the result in same unit
        if self.unit == other.unit:
            return ((self.value - other.value), self.unit)
        # Otherwise, return result per second
        values = self._standardize_units(other, target_unit='second')
        diff_per_second = values[0] - values[1]
        # result is per second
        return (diff_per_second, 'second')

    def __gt__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value > other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] > values[1])

    def __ge__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value >= other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] >= values[1])

    def __lt__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value < other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] < values[1])

    def __le__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value <= other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] <= values[1])

    def __eq__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value == other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] == values[1])

    def __ne__(self, other):
        # Only proceed if entities are the same
        if not self.entity == other.entity:
            raise ValueError("Entities have to be the same for comparison")
        # Compare values right away if units are the same
        if self.unit == other.unit:
            return self.value != other.value
        # If units are different, convert to seconds first
        values = self._standardize_units(other, target_unit='second')
        return (values[0] != values[1])

    # Unit conversions
    @staticmethod
    def convert_unit(source_unit, target_unit, source_value,
                     source_period=None, target_period=None):
        """Convert value per unit from source to target unit. If a unit is
        absolute, total timedelta period has to be provided. If a unit is a
        month or a year, it is recommended to pass timedelta period object
        directly, if not provided, the approximation will be used.
        """
        if (target_unit == 'absolute' and not target_period) or \
           (source_unit == 'absolute' and not source_period):
                raise ValueError('Absolute period needs to be provided.')
                return
        if not source_period:
            source_period = QuantitativeState._get_period_from_unit(
                source_unit)
        if not target_period:
            target_period = QuantitativeState._get_period_from_unit(
                target_unit)
        if source_unit == 'second':
            return QuantitativeState.from_seconds(source_value, target_period)
        if target_unit == 'second':
            return QuantitativeState.value_per_second(
                source_value, source_period)
        return round(source_value * (target_period / source_period))

    @staticmethod
    def value_per_second(value, period):
        """Get value per second given total value per period and a timedelta
        period object."""
        seconds = period.total_seconds()
        per_second = value / seconds
        return per_second

    @staticmethod
    def from_seconds(value_per_second, period):
        """Get total value per given period given timedelta period object and
        value per second."""
        seconds = period.total_seconds()
        total_value = value_per_second * seconds
        return total_value

    def _standardize_units(self, other, target_unit='second',
                           target_period=None):
        values = []
        for state in [self, other]:
            if state.unit != target_unit:
                new_value = self.convert_unit(source_unit=state.unit,
                                              target_unit=target_unit,
                                              source_value=state.value,
                                              target_period=target_period)
            values.append(new_value)
        return values

    @staticmethod
    def _get_period_from_unit(unit):
        if unit == 'day':
            return timedelta(days=1)
        if unit == 'week':
            return timedelta(weeks=1)
        if unit == 'hour':
            return timedelta(hours=1)
        if unit == 'minute':
            return timedelta(minutes=1)
        if unit == 'second':
            return timedelta(seconds=1)
        if unit == 'month':
            logger.warning('Using approximate value 30 days for a month.')
            return timedelta(days=30)
        if unit == 'year':
            logger.warning('Using approximate value 365 days for a year.')
            return timedelta(days=365)
        if unit == 'absolute':
            logger.warning('Cannot derive absolute period, '
                           'it needs to be provided.')
            return None
