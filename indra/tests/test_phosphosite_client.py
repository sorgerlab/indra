from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases.phosphosite_client import map_to_human_site
from nose.plugins.attrib import attr


def test_map_mouse_to_human():
    mouse_up_id = 'Q61337'
    pos = map_to_human_site(mouse_up_id, 'S', '112')
    assert pos == '75'


def test_isoform_mapping_from_human():
    up_id = 'P29353'
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


def test_isoform_mapping_from_mouse():
    up_id = 'P29353'
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


if __name__ == '__main__':
    test_isoform_mapping_from_human()
