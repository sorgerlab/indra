from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.sources.lincs.api import _get_lincs_drug_target_data, \
    process_from_web


def test_get_drug_target_data():
    data_list = _get_lincs_drug_target_data()
    assert len(data_list) > 100, len(data_list)


def test_process_from_web():
    lincs_p = process_from_web()
    assert lincs_p is not None
    assert lincs_p.statements
