from indra.reach import reach_api
from indra.reach.processor import ReachProcessor

def test_parse_site_text():
    rp = ReachProcessor(None)
    text = ['threonine 185', 'thr 185', 'thr-185', 'threonine residue 185', 'T185']
    for t in text:
        residue, site = rp._parse_site_text(t)
        assert(residue == 'Threonine')
        assert(site == '185')

def test_parse_site_residue_only():
    rp = ReachProcessor(None)
    text = ['serine residue', 'serine', 'a serine site']
    for t in text:
        residue, site = rp._parse_site_text(t)
        assert(residue == 'Serine')
        assert(site is None)

def test_read_sentence_online():
    rp = reach_api.process_text('MEK1 phosphorylates ERK2.')
    assert(len(rp.statements) == 1)
