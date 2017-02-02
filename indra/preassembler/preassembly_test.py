if __name__ == '__main__':
    from indra.tools import assemble_corpus as ac
    import sys

    stmts = ac.load_statements(sys.argv[1])
    rel = ac.run_preassembly(stmts)

