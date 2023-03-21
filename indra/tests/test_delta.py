from indra.statements.delta import *
import pytest
from datetime import date, timedelta


def test_qualitative_serialization():
    qd = QualitativeDelta(polarity=1, adjectives=['strong', 'significant'])
    qdj = qd.to_json()
    assert qdj['type'] == 'qualitative'
    assert qdj['polarity'] == 1
    assert qdj['adjectives'] == ['strong', 'significant']
    assert qdj == Delta.from_json(qdj).to_json()


def test_change_parameters():
    qd = QualitativeDelta(polarity=None, adjectives=None)
    assert qd.polarity is None
    assert qd.adjectives == []
    # Set polarity
    qd.set_polarity(1)
    assert qd.polarity == 1
    qd.set_polarity(-1)
    assert qd.polarity == -1
    qd.set_polarity(None)
    assert qd.polarity is None
    # Add adjectives
    qd.add_adjectives(['strong'])
    assert qd.adjectives == ['strong']
    qd.add_adjectives(['significant', 'severe'])
    assert qd.adjectives == ['strong', 'significant', 'severe']


def test_qualitative_delta_comparisons():
    qd1 = QualitativeDelta(polarity=1, adjectives=['strong', 'significant'])
    qd2 = QualitativeDelta(polarity=1, adjectives=['strong', 'significant'])
    qd3 = QualitativeDelta(polarity=-1, adjectives=['strong', 'significant'])
    qd4 = QualitativeDelta(polarity=1, adjectives=['weak', 'insignificant'])
    qd5 = QualitativeDelta(polarity=None, adjectives=None)
    assert qd1.equals(qd2)
    assert not qd1.equals(qd3)
    assert not qd1.equals(qd4)
    assert qd1.is_opposite(qd3)
    assert not qd1.is_opposite(qd5)
    assert not qd1.is_opposite(qd4)


def test_quantitative_serialization():
    qs = QuantitativeState('person', 100, 'absolute', 'more than',
                           'More than 100 people arrived to Ethiopia.', 1)
    qsj = qs.to_json()
    assert qsj['type'] == 'quantitative'
    assert qsj['entity'] == 'person', qsj['entity']
    assert qsj['value'] == 100, qsj['value']
    assert qsj['unit'] == 'absolute', qsj['unit']
    assert qsj['modifier'] == 'more than', qsj['modifier']
    assert qsj['polarity'] == 1
    assert qsj == Delta.from_json(qsj).to_json()


def test_quantitative_conversions():
    day = QuantitativeState._get_period_from_unit('day')
    assert QuantitativeState.value_per_second(86400, day) == 1
    assert QuantitativeState.from_seconds(1, day) == 86400
    qs1 = QuantitativeState(value=86400, unit='day')
    qs2 = QuantitativeState(value=302400, unit='week')
    values = qs1._standardize_units(qs2, target_unit='second')
    assert values[0] == 1.0  # 86400 / 86400 sec/day
    assert values[1] == 0.5  # 302400 / 604800 sec/week
    # Convert between different rates
    assert QuantitativeState.convert_unit('day', 'week', 10) == 70
    assert QuantitativeState.convert_unit('week', 'day', 21) == 3
    # Use given periods versus approximate values
    jan = date(2019, 2, 1) - date(2019, 1, 1)
    feb = date(2019, 3, 1) - date(2019, 2, 1)
    feb_leap = date(2016, 3, 1) - date(2016, 2, 1)
    assert QuantitativeState.convert_unit(
        'day', 'month', 1, target_period=jan) == 31
    assert QuantitativeState.convert_unit(
        'day', 'month', 1, target_period=feb) == 28
    assert QuantitativeState.convert_unit(
        'day', 'month', 1, target_period=feb_leap) == 29
    assert QuantitativeState.convert_unit('day', 'month', 1) == 30
    # Convert to absolute value
    abs_period = timedelta(weeks=5)
    assert QuantitativeState.convert_unit(
        'day', 'absolute', 2, day, abs_period) == 70
    assert QuantitativeState.convert_unit(
        'week', 'absolute', 25, target_period=abs_period) == 125
    # Get rate from absolute value
    assert QuantitativeState.convert_unit(
        'absolute', 'day', 70, source_period=abs_period) == 2
    assert QuantitativeState.convert_unit(
        'absolute', 'week', 70, source_period=abs_period) == 14
    # Convert to or from absolute value without providing a total period
    with pytest.raises(ValueError):
        QuantitativeState.convert_unit('absolute', 'week', 70)
    with pytest.raises(ValueError):
        QuantitativeState.convert_unit('day', 'absolute', 2)


def test_arithmetic_operations():
    qs1 = QuantitativeState('person', 15, 'day')
    qs2 = QuantitativeState('person', 10, 'day')
    qs3 = QuantitativeState('person', 100, 'week')
    qs4 = QuantitativeState('person', 70, 'week')
    qs5 = QuantitativeState('box', 100, 'day')
    qs6 = QuantitativeState('person', 100, 'absolute')
    qs7 = QuantitativeState('person', 30, 'absolute')
    day = QuantitativeState._get_period_from_unit('day')
    # Operations with the same entity and unit
    assert qs1 + qs2 == (25, 'day'), qs1 + qs2
    assert qs1 - qs2 == (5, 'day'), qs1 - qs2
    assert qs1 > qs2
    assert qs2 < qs1
    assert qs1 != qs2
    # Operations with absolute values
    assert qs6 + qs7 == (130, 'absolute')
    assert qs6 - qs7 == (70, 'absolute')
    assert qs6 > qs7
    assert qs7 < qs6
    # Operations with different units
    sum_per_second = (qs1 + qs4)[0]
    sum_per_day = QuantitativeState.from_seconds(sum_per_second, day)
    assert sum_per_day == 25, sum_per_day
    diff_per_second = (qs1 - qs4)[0]
    diff_per_day = QuantitativeState.from_seconds(diff_per_second, day)
    assert diff_per_day == 5, diff_per_day
    assert qs1 > qs4
    assert qs2 < qs3
    assert qs2 == qs4
    assert qs1 != qs3
    # Operations with different entities
    with pytest.raises(ValueError):
        qs1 + qs5
    # Operations with absolute unit when absolute period not provided
    with pytest.raises(ValueError):
        qs1 + qs6
