from indra.reach import reach_api
from indra.reach.processor import ReachProcessor

def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'Threonine')
        assert(site == '185')

def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site']
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'Serine')
        assert(site is None)

def test_phosphorylate():
    rp = reach_api.process_text('MEK1 phosphorylates ERK2.')
    assert(len(rp.statements) == 1)

def test_bind():
    rp = reach_api.process_text('MEK1 binds ERK2.')
    assert(len(rp.statements) == 1)

def test_activity():
    rp = reach_api.process_text('MEK1 activates ERK2.')
    assert(len(rp.statements) == 1)
