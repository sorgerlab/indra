from indra.statements.delta import *
from nose.tools import assert_raises


def test_quantitative_serialization():
    qs = QuantitativeState('person', 100, 'absolute', 'more than',
                           'More than 100 people arrived to Ethiopia.')
    qsj = qs.to_json()
    assert qsj['entity'] == 'person', qsj['entity']
    assert qsj['value'] == 100, qsj['value']
    assert qsj['unit'] == 'absolute', qsj['unit']
    assert qsj['modifier'] == 'more than', qsj['modifier']
    assert qsj == Delta.from_json(qsj).to_json()


def test_quantitative_conversions():
    qs = QuantitativeState()
    assert qs.value_per_second(86400, 1) == 1
    assert qs.from_seconds(1, 1) == 86400
    qs2 = QuantitativeState(value=86400, unit='daily')
    qs3 = QuantitativeState(value=302400, unit='weekly')
    values = qs2._standardize_units(qs3, target_unit='per_second')
    assert values[0] == 1.0  # 86400 / 86400 sec/day
    assert values[1] == 0.5  # 302400 / 604800 sec/week


def test_arithmetic_operations():
    qs1 = QuantitativeState('person', 15, 'daily')
    qs2 = QuantitativeState('person', 10, 'daily')
    qs3 = QuantitativeState('person', 100, 'weekly')
    qs4 = QuantitativeState('person', 70, 'weekly')
    qs5 = QuantitativeState('box', 100, 'daily')
    # Operations with the same entity and unit
    assert (qs1 + qs2) == 25
    assert (qs1 - qs2) == 5
    assert qs1 > qs2
    assert qs2 < qs1
    assert qs1 != qs2
    # Operations with different units
    assert (qs1 + qs4) == 25  # daily
    assert (qs1 - qs4) == 5  # daily
    assert qs1 > qs4
    assert qs2 < qs3
    assert qs2 == qs4
    assert qs1 != qs3
    # Operations with different entities
    assert_raises(ValueError, qs1.__add__, qs5)
