import pickle
import numpy
import time
from pysb.integrate import Solver
from belpy.pysb_assembler import *


if __name__ == '__main__':
    # Instantiate the assembler
    pa = PysbAssembler()
    
    # TRIPS processing   
    use_xml = False
    if use_xml:
        fname = 'test.xml'
        print 'Processing parser output from XML file %s...' % fname
        tp = trips_api.process_xml(open(fname).read())
    else:
        print 'Submitting text to TRIPS parser...'
        tstart = time.time()
        print 'Processing parser output...' % fname
        tp = trips_api.process_text(open('model_text.txt').read())
        tend = time.time()
        print 'TRIPS parser took %d seconds.' % (tend - tstart)
    pa.add_statements(tp.belpy_stmts)

    # BioPAX processing
    use_owl = True
    if use_owl:
        fname = 'DUSP.owl'
        print 'Processing OWL file %s' % fname
        bp = biopax_api.process_owl('DUSP.owl')
    else:
        print 'Processing output from PathwayCommons query'
        bp = biopax_api.process_pc_neighborhood(['MAPK1'])
    bp.get_dephosphorylation(force_contains=['MAPK1'])
    pa.add_statements(bp.belpy_stmts)
    
    # Assemble model
    model = pa.make_model()

    #ts = numpy.linspace(0,10,10)
    #solver = Solver(model, ts)
    #solver.run()
