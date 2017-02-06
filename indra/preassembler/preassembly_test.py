if __name__ == '__main__':
    from indra.tools import assemble_corpus as ac
    import logging
    import sys

    aclogger = logging.getLogger('assemble_corpus')
    aclogger.setLevel(logging.DEBUG)

    prelogger = logging.getLogger('preassembler')
    prelogger.setLevel(logging.DEBUG)

    stmts = ac.load_statements(sys.argv[1])
    rel = ac.run_preassembly(stmts, poolsize=None)

