from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.databases import biogrid_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr
from indra.sources.biogrid import process_file
from indra.statements import Complex


@attr('webservice', 'nonpublic')
def test_biogrid_request():
    results = biogrid_client._send_request(['MAP2K1', 'MAPK1'])
    assert results is not None
    assert unicode_strs(results)


def test_biogrid_tsv():
    # Download biogrid file form the web and process it
    bp = process_file(None)

    # We should have a lot of statementse
    statements = bp.statements
    assert(len(statements) > 500000)

    # Any given statement should be a complex, with appropriate evidence
    s0 = statements[0]
    assert(isinstance(s0, Complex))
    ev = s0.evidence[0]
    assert(ev.source_api == 'biogrid')
    assert(ev.text is None)
    assert(ev.pmid is not None)
    assert('tsv_row' in ev.annotations)

    # The first statement in the file involves MAP2K4 and FLNC
    assert(str(s0.members[0]) == 'MAP2K4()')
    assert(str(s0.members[1]) == 'FLNC()')
