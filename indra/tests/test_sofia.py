from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
from nose.plugins.attrib import attr

from indra.sources import sofia


@attr('webservice', 'nonpublic')
def test_text_process_webservice():
    txt = 'rainfall causes floods'
    sp = sofia.process_text(txt)
    assert len(sp.statements) == 1
    assert sp.statements[0].subj.name == 'rainfall'
    assert sp.statements[0].obj.name == 'floods'


def test_process_json():
    path_here = os.path.abspath(os.path.dirname(__file__))
    test_file = os.path.join(path_here, 'sofia_test.json')
    with open(test_file, 'r') as fh:
        js = json.load(fh)
    sp = sofia.process_json(js)
    assert len(sp.statements) == 1
    assert sp.statements[0].subj.name == 'rainfall'
    assert sp.statements[0].obj.name == 'floods'
