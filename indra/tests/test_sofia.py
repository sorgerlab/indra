from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import os
import json
from nose.plugins.attrib import attr

from indra.sources import sofia
from indra.statements.statements import Influence, Event
from indra.statements.context import WorldContext


# Tell nose to not run tests in the imported modules
Influence.__test__ = False
Event.__test__ = False
WorldContext.__test__ = False


@attr('webservice', 'nonpublic')
def test_text_process_webservice():
    txt = 'rainfall causes floods'
    sp = sofia.process_text(txt)
    assert len(sp.statements) == 1
    assert sp.statements[0].subj.concept.name == 'rainfall'
    assert sp.statements[0].obj.concept.name == 'floods'


def test_process_json():
    path_here = os.path.abspath(os.path.dirname(__file__))
    test_file = os.path.join(path_here, 'sofia_test.json')
    with open(test_file, 'r') as fh:
        js = json.load(fh)
    sp = sofia.process_json(js)
    assert len(sp.statements) == 2
    assert isinstance(sp.statements[0], Influence)
    assert sp.statements[0].subj.concept.name == 'rainfall'
    assert sp.statements[0].obj.concept.name == 'floods'
    assert len(sp.statements[0].evidence) == 1, len(sp.statements[0].evidence)
    assert isinstance(sp.statements[1], Event)
    assert sp.statements[1].concept.name == 'inflation'
    assert isinstance(sp.statements[1].context, WorldContext)
    assert sp.statements[1].context.time.text == '28, JULY, 2016'
    assert sp.statements[1].context.geo_location.name == 'South Sudan'
