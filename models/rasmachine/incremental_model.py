import pickle
from indra.assemblers import PysbAssembler

class IncrementalModel(object):
    def __init__(self, fname=None):
        if fname is None:
            self.stmts = {}
        else:
            try:
                self.stmts = pickle.load(open(fname, 'rb'))
            except:
                print 'Could not load %s' % fname
                self.stmts = {}

    def save(self, fname='model.pkl'):
        with open(fname, 'wb') as fh:
            pickle.dump(self.stmts, fh)

    def add_statements(self, pmid, stmts):
        self.stmts[pmid] = stmts

    def get_statements(self):
        stmt_lists = [v for k, v in self.stmts.iteritems()]
        stmts = []
        for s in stmt_lists:
            stmts += s
        return stmts

    def make_model(self):
        pa = PysbAssembler()
        pa.add_statements(self.get_statements())
        pa.make_model()
        return pa.model
