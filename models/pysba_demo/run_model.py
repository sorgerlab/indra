import pickle
import numpy
import time
from pysb.integrate import Solver
from belpy.pysb_assembler import *


if __name__ == '__main__':
    pa = PysbAssembler()
    use_pickle = False
    if use_pickle:
        tp = pickle.load(open('model_text_tp.pkl'))
    else:
        tstart = time.time()
        tp = trips_api.process_text(open('model_text.txt').read())
        tend = time.time()
        print 'TRIPS parser took %d seconds.' % (tend - tstart)
    bp = biopax_api.process_owl('DUSP.owl')
    bp.get_dephosphorylation()
    pa.add_statements(tp.belpy_stmts)
    pa.add_statements(bp.belpy_stmts)
    model = pa.make_model()
    
    ts = numpy.linspace(0,10,10)
    solver = Solver(model, ts)
    solver.run()
