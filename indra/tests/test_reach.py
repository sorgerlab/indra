from indra import reach
from indra.reach.processor import ReachProcessor

def test_parse_site_text():
    text = ['threonine 185', 'thr 185', 'thr-185',
            'threonine residue 185', 'T185']
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'T')
        assert(site == '185')

def test_parse_site_residue_only():
    text = ['serine residue', 'serine', 'a serine site']
    for t in text:
        residue, site = ReachProcessor._parse_site_text(t)
        assert(residue == 'S')
        assert(site is None)

def test_valid_name():
    assert(ReachProcessor._get_valid_name('') == '')
    assert(ReachProcessor._get_valid_name('a') == 'a')
    assert(ReachProcessor._get_valid_name('Name123') == 'Name123')
    assert(ReachProcessor._get_valid_name('<>#~!,./][;-') == '____________')
    assert(ReachProcessor._get_valid_name('PI3 Kinase') == 'PI3_Kinase')
    assert(ReachProcessor._get_valid_name('14-3-3') == 'p14_3_3')

def test_phosphorylate():
    rp = reach.process_text('MEK1 phosphorylates ERK2.')
    assert(len(rp.statements) == 1)
    s = rp.statements[0]
    assert (s.enz.name == 'MAP2K1')
    assert (s.sub.name == 'MAPK1')

def test_activate():
    rp = reach.process_text('HRAS activates BRAF.')
    assert(len(rp.statements) == 1)
    s = rp.statements[0]
    assert (s.subj.name == 'HRAS')
    assert (s.obj.name == 'BRAF')

def test_bind():
    rp = reach.process_text('MEK1 binds ERK2.')
    assert(len(rp.statements) == 1)

def test_activity():
    rp = reach.process_text('MEK1 activates ERK2.')
    assert(len(rp.statements) == 1)

def test_mutation():
    rp = reach.process_text('BRAF(V600E) phosphorylates MEK.')
    assert(len(rp.statements) == 1)
    braf = rp.statements[0].enz
    assert(braf.name == 'BRAF')
    assert(len(braf.mutations) == 1)
    assert(braf.mutations[0].position == '600')
    assert(braf.mutations[0].residue_from == 'V')
    assert(braf.mutations[0].residue_to == 'E')
