from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases.phosphosite_client import map_to_human_site

def test_map_mouse_to_human():
    mouse_up_id = 'Q61337'
    (res, pos) = map_to_human_site(mouse_up_id, 'S', '112')
    assert res == 'S'
    assert pos == '75'

if __name__ == '__main__':
    test_map_mouse_to_human()
