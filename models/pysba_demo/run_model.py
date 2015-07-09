import pickle
import numpy
import time
import argparse
from pysb.integrate import Solver
import pysb.tools.render_reactions as rr
import pysb.tools.render_species as rs
from belpy.pysb_assembler import *


if __name__ == '__main__':
    # Instantiate the assembler
    pa = PysbAssembler()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--offline', action="store_true")
    parser.add_argument('--render', action="store_true")
    args = parser.parse_args()

    if args.offline:
        use_xml = True
        use_owl = True
    else:
        use_xml = False
        use_owl = False


    # TRIPS processing   
    if use_xml:
        fname = 'test.xml'
        print 'Processing parser output from XML file %s...' % fname
        tp = trips_api.process_xml(open(fname).read())
    else:
        print 'Submitting text to TRIPS parser...'
        tstart = time.time()
        tp = trips_api.process_text(open('model_text.txt').read())
        tend = time.time()
        print 'TRIPS parser took %d seconds.' % (tend - tstart)
    pa.add_statements(tp.belpy_stmts)

    # BioPAX processing
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

    # Model updates
    Egfr = model.monomers['Egfr']
    r1 = model.rules[u'Egfr_phospho_Egfr_Y1068']
    r1.rule_expression = (Egfr(Egfr=1) % Egfr(Egfr=1, Y1068='u') >> Egfr(Egfr=1) % Egfr(Egfr=1, Y1068='p'))
    r2 = model.rules[u'Egfr_phospho_Egfr_Y1148']
    r2.rule_expression = (Egfr(Egfr=1) % Egfr(Egfr=1, Y1148='u') >> Egfr(Egfr=1) % Egfr(Egfr=1, Y1148='p'))
    for r in (r1, r2):
        r.reactant_pattern = r.rule_expression.reactant_pattern
        r.product_pattern = r.rule_expression.product_pattern

    if args.render:
        print 'Rendering reactions...'
        gv_str = rr.run(model)
        with open('reactions.dot', 'w') as f:
            f.write(gv_str)

        print 'Rendering species...'
        gv_str = rs.run(model)
        with open('species.dot', 'w') as f:
            f.write(gv_str)

      

    # Run model simulation
    # TODO: save figure as file
    #ts = numpy.linspace(0,10,10)
    #solver = Solver(model, ts)
    #solver.run()
