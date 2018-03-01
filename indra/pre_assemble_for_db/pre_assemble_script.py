import indra.tools.assemble_corpus as ac

def process_statements(stmts):
    stmts = ac.map_grounding(stmts)
    stmts = ac.map_sequence(stmts)
    stmts = ac.run_preassembly(stmts, return_toplevel=False)
    return stmts
