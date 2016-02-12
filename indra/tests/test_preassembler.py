from indra.preassembler import Preassembler
from indra.trips import trips_api
from indra.statements import Agent, Phosphorylation, BoundCondition

def test_from_text():
    sentences = ['Src phosphorylates Ras, bound to GTP, at Tyr32.', 
                'Src phosphorylates NRAS at Tyr32.',
                'Src phosphorylates NRAS that is bound to GTP.']
    for s in sentences:
        tp = trips_api.process_text(s)
        pa.add_statements(tp.statements)

def test_duplciates():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    st1 = Phosphorylation(src, ras, 'Phosphorylation', None)
    st2 = Phosphorylation(src, ras, 'Phosphorylation', None)
    pa = Preassembler([st1, st2])
    stmts = pa.assemble()
    assert(len(stmts) == 1)

def test_src_phos_nras():
    src = Agent('SRC', db_refs = {'HGNC': '11283'})
    ras = Agent('RAS', db_refs = {'FA': '03663'})
    gtp = Agent('GTP', db_refs = {'CHEBI': '15996'})
    rasgtp = Agent('RAS', db_refs = {'FA': '03663'}, 
        bound_conditions=[BoundCondition(gtp, True)])
    nras = Agent('NRAS', db_refs = {'HGNC': '7989'})
    nrasgtp = Agent('NRAS', db_refs = {'HGNC': '7989'}, 
        bound_conditions=[BoundCondition(gtp, True)])
    st1 = Phosphorylation(src, rasgtp, 'PhosphorylationTyrosine', '32')
    st2 = Phosphorylation(src, nras, 'PhosphorylationTyrosine', '32')
    st3 = Phosphorylation(src, nrasgtp, 'Phosphorylation', None)
    pa = Preassembler([st1, st2, st3])
    stmts = pa.assemble()

    assert(len(stmts) == 1)
    assert(stmts[0].enz.name == 'Src')
    assert(stmts[0].sub.name == 'NRAS')
    assert(stmts[0].mod = 'PhosphorylationTyrosine')
    assert(stmts[0].mod_pos = '32')
