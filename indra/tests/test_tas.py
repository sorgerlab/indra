from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.tas.api import _load_data


def test_load_data():
    data = _load_data()
    assert len(data) > 100, len(data)
