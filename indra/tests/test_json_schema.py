import json
import jsonschema
import os
from nose.tools import assert_raises
from jsonschema.exceptions import ValidationError

dir_this = os.path.dirname(__file__)
schema_file = os.path.join(dir_this, '../../schemas/statements_schema.json')
with open(schema_file, 'r') as f:
    schema = json.loads(f.read())

valid_agent1 = {'name': 'RAF', 'db_refs': {'TEXT': 'RAF'}}
valid_agent2 = {'name': 'RAS', 'db_refs': {'TEXT': 'RAS'}}

invalid_agent1 = {'name': 'RAS', 'db_refs': 2}
invalid_agent2 = {'db_refs': {'TEXT': 'RAS'}}
invalid_agent3 = {'name': 'cow', 'db_refs': {'TEXT': 'RAS'},
                  'bound_conditions': 'mooooooooooo!'}


def val(s):
    jsonschema.validate([s], schema)


def test_valid_phosphorylation():
    # Calls to validate() in this function should not raise exceptions
    s = {'enz': valid_agent1, 'sub': valid_agent2, 'type': 'Phosphorylation',
         'id': '5'}
    jsonschema.validate([s], schema)

    s['residue'] = 'S'
    jsonschema.validate([s], schema)

    s['position'] = '10'
    jsonschema.validate([s], schema)


def test_invalid_phosphorylation():
    s = {'enz': valid_agent1, 'sub': valid_agent2, 'type': 'Phosphorylation',
         'id': '5', 'residue': 5}  # residue should be a string
    assert_raises(ValidationError, val, s)

    s = {'enz': valid_agent1, 'sub': invalid_agent1, 'type': 'Phosphorylation',
         'id': '5'}
    assert_raises(ValidationError, val, s)

    s = {'enz': valid_agent1, 'sub': invalid_agent2, 'type': 'Phosphorylation',
         'id': '5'}
    assert_raises(ValidationError, val, s)

    s = {'enz': valid_agent1, 'sub': invalid_agent3, 'type': 'Phosphorylation',
         'id': '5'}
    assert_raises(ValidationError, val, s)


def test_valid_active_form():
    s = {'agent': valid_agent1, 'activity': 'kinase', 'is_active': True,
         'type': 'ActiveForm', 'id': '6'}
    jsonschema.validate([s], schema)


def test_invalid_active_form():
    s = {'agent': invalid_agent1, 'activity': 'kinase', 'is_active': True,
         'type': 'ActiveForm', 'id': '6'}
    assert_raises(ValidationError, val, s)

    s = {'agent': valid_agent1, 'activity': 'kinase', 'is_active': 'moo',
         'type': 'ActiveForm', 'id': '6'}
    assert_raises(ValidationError, val, s)

    s = {'agent': valid_agent1, 'activity': 'kinase', 'is_active': True,
         'type': 'ActiveForm', 'id': 42}
    assert_raises(ValidationError, val, s)

    s = {'agent': valid_agent1, 'activity': 'kinase', 'is_active': True,
         'type': 'MOO', 'id': '6'}
    assert_raises(ValidationError, val, s)

    s = {'agent': valid_agent1, 'activity': {'cow': False}, 'is_active': True,
         'type': 'ActiveForm', 'id': '6'}
    assert_raises(ValidationError, val, s)


def test_valid_complex():
    s = {'members': [valid_agent1, valid_agent2], 'type': 'Complex', 'id': '3'}
    jsonschema.validate([s], schema)

    s = {'members': [], 'type': 'Complex', 'id': '3'}
    jsonschema.validate([s], schema)


def test_invalid_complex():
    s = {'members': [invalid_agent1, valid_agent2], 'type': 'Complex',
         'id': '3'}
    assert_raises(ValidationError, val, s)

    s = {'members': [valid_agent1, invalid_agent2], 'type': 'Complex',
         'id': '3'}
    assert_raises(ValidationError, val, s)
