from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
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
    sp = sofia.process_json()
    assert len(sp.statements) == 1
    assert sp.statements[0].subj.name == 'rainfall'
    assert sp.statements[0].obj.name == 'floods'
