from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from nose.plugins.attrib import attr

from indra.sources import sofia


@attr('webservice', 'nonpublic')
def test_text_process_webservice():
    txt = 'rainfall causes floods'
    sp = sofia.process_text(txt)
    assert len(sp.statements) == 1
    assert sp.statements[0].subj.name == 'rainfall'
    assert sp.statements[0].obj.name == 'floods'
