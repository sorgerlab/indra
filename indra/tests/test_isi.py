from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.sources.isi.api import process_text
from indra.statements import Complex, Phosphorylation


def test_process_complex():
    ip = process_text('Ras binds to Raf.', '42')

    statements = ip.statements
    assert len(statements) == 1

    s0 = statements[0]
    assert type(s0) == Complex

    assert len(s0.members) == 2

    m0 = s0.members[0]
    assert m0.name == 'RAS'
    assert m0.db_refs == {'TEXT': 'Ras', 'PF': 'PF00071', 'FPLX': 'RAS'}

    m1 = s0.members[1]
    assert m1.name == 'RAF'
    assert m1.db_refs == {'TEXT': 'Raf', 'FPLX': 'RAF'}

    assert len(s0.evidence) == 1
    ev = s0.evidence[0]
    assert ev.source_api == 'isi'
    assert ev.pmid == '42'
    assert ev.text == 'Ras binds to Raf.'
    assert ev.annotations['interaction'] == ['binds', None, 'Ras', 'Raf']
    assert ev.annotations['source_id'] is not None


def test_process_phosphorylation():
    # Include a sentence without a mechanism to test ISI's ability to
    # associate the relevant sentence text with an extracted event.
    ip = process_text('This sentence is false. Ras phosphorylates Raf.', '42')

    statements = ip.statements
    assert len(statements) == 1

    s0 = statements[0]
    assert type(s0) == Phosphorylation

    enz = s0.enz
    assert enz.name == 'RAS'
    assert enz.db_refs == {'TEXT': 'Ras', 'PF': 'PF00071', 'FPLX': 'RAS'}

    sub = s0.sub
    assert sub.name == 'RAF'
    assert sub.db_refs == {'TEXT': 'Raf', 'FPLX': 'RAF'}

    assert len(s0.evidence) == 1
    ev = s0.evidence[0]
    assert ev.source_api == 'isi'
    assert ev.source_id is not None
    assert ev.pmid == '42'
    assert ev.text == 'Ras phosphorylates Raf.'
    assert ev.annotations['interaction'] == ['phosphorylates', 'Ras', 'Raf']
    assert ev.annotations['source_id'] is not None
