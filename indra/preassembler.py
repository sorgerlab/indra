
class Preassembler(object):

    def __init__(self, stmts=None):
        if stmts:
            self.stmts = stmts
        else:
            self.stmts = []

    def add_statements(self, stmts):
        """Add to the current list of statements."""
        self.stmts += stmts

    def assemble(self):
        for s in self.stmts:
            print s

    def get_statement_family(stmt):
        pass

if __name__ == '__main__':
    import pickle
    with open('example_kami.pck') as f:
        stmts = pickle.load(f)

    pa = Preassembler(stmts)
    pa.assemble()
