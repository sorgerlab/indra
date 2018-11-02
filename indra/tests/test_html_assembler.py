from indra.statements import *
from indra.assemblers.html import HtmlAssembler

def make_stmt():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FPLX': 'RAS'})
    ev = Evidence(text="We noticed that the Src kinase was able to "
                       "phosphorylate Ras proteins.",
                  source_api='test', pmid='1234567',
                  annotations={'agents': {'raw_text': ['Src kinase',
                                                       'Ras proteins']}})
    st = Phosphorylation(src, ras, 'tyrosine', '32', evidence=[ev])
    return st

def test_format_evidence_text():
    stmt = make_stmt()
    ev_list = HtmlAssembler.format_evidence_text(stmt)
    assert len(ev_list) == 1
    ev = ev_list[0]
    assert isinstance(ev, dict)
    assert set(ev.keys()) == set(['source_api', 'pmid', 'text'])
    assert ev['source_api'] == 'test'
    assert ev['pmid'] == '1234567'
    assert ev['text'] == ('We noticed that the '
                          '<span class="label label-subject">Src kinase</span> '
                          'was able to phosphorylate '
                          '<span class="label label-object">'
                          'Ras proteins</span>.')
