import pickle
import numpy
import time
import argparse
import matplotlib.pyplot as plt
from pysb.integrate import Solver
import pysb.tools.render_reactions as rr
import pysb.tools.render_species as rs
from indra import trips, biopax
from indra.assemblers import PysbAssembler


def set_initial(model, cell_line):
    try:
        x0 = numpy.recfromcsv('total_prots.csv') 
        x0_names = x0['proteins']
        x0_values = x0[cell_line.lower()]
        print 'Setting initial conditions for cell line %s...' % cell_line
        for n, v in zip(x0_names, x0_values):
            try:
                monomer = model.monomers[n]
            except KeyError:
                continue
            set_base_initial_condition(model, monomer, v)
    except IOError:
        print 'Could not open total proteins file.'

if __name__ == '__main__':
    # Instantiate the assembler
    pa = PysbAssembler()

    parser = argparse.ArgumentParser()
    parser.add_argument('--online', action="store_true")
    parser.add_argument('--render', action="store_true")
    parser.add_argument('--simulate', action="store_true")
    parser.add_argument('--loadmodel', action="store_true")
    args = parser.parse_args()

    if args.online:
        use_xml = False
        use_owl = False
    else:
        use_xml = True
        use_owl = True


    # TRIPS processing   
    if use_xml:
        fname = 'egfr_model_v4.xml'
        print 'Processing parser output from XML file %s...' % fname
        tp = trips.process_xml(open(fname).read())
    else:
        print 'Submitting text to TRIPS parser...'
        tstart = time.time()
        tp = trips.process_text(open('model_text.txt').read())
        tend = time.time()
        print '> TRIPS parser took %d seconds.' % (tend - tstart)
    pa.add_statements(tp.statements)

    # BioPAX processing
    if use_owl:
        fname = 'DUSP.owl'
        print 'Processing OWL file %s...' % fname
        bp = biopax.process_owl('DUSP.owl')
    else:
        print 'Processing output from PathwayCommons query'
        bp = biopax.process_pc_neighborhood(['DUSP4'])
    bp.get_dephosphorylation(force_contains=['MAPK1'])
    pa.add_statements(bp.statements)

    # Assemble model
    model = pa.make_model(initial_conditions=True)


    RAF1 = model.monomers['RAF1']
    HRAS = model.monomers['HRAS']

    r1 = model.rules['RAF1_autophospho_RAF1_T491']
    r1.rule_expression = (RAF1(HRAS=1, T491='u') % HRAS(RAF1=1) >> RAF1(HRAS=1, T491='p') % HRAS(RAF1=1))
    r2 = model.rules['RAF1_autophospho_RAF1_S494']
    r2.rule_expression = (RAF1(HRAS=1, S494='u') % HRAS(RAF1=1) >> RAF1(HRAS=1, S494='p') % HRAS(RAF1=1))
    for r in (r1, r2):
        r.reactant_pattern = r.rule_expression.reactant_pattern
        r.product_pattern = r.rule_expression.product_pattern


    # Rendering reactions and species
    if args.render:
        print 'Rendering reactions...'
        gv_str = rr.run(model)
        with open('reactions.dot', 'w') as f:
            f.write(gv_str)

        print 'Rendering species...'
        gv_str = rs.run(model)
        with open('species.dot', 'w') as f:
            f.write(gv_str)

    if args.simulate:
        # Set parameters
        print 'Setting model parameters...'
        model.add_component(Parameter("kp1", 1.667e-06))  # ligand-monomer binding (scaled), units: /molecule/s
        model.add_component(Parameter("km1", 0.06))       # ligand-monomer dissociation, units: /s
        model.add_component(Parameter("kp2", 5.556e-06))  # aggregation of bound monomers (scaled), units: /molecule/s
        model.add_component(Parameter("km2", 0.1))        # dissociation of bound monomers, units: /s
        model.add_component(Parameter("kp3", 0.5))        # dimer transphosphorylation, units: /s
        model.add_component(Parameter("km3", 4.505))      # dimer dephosphorylation, units: /s
        model.add_component(Parameter("kp14", 3))         # Shc transphosphorylation, units: /s
        model.add_component(Parameter("km14", 0.03))      # Shc dephosphorylation, units: /s
        model.add_component(Parameter("km16", 0.005))     # Shc cytosolic dephosphorylation, units: /s
        model.add_component(Parameter("kp9", 8.333e-07))  # binding of Grb2 to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km9", 0.05))       # dissociation of Grb2 from receptor, units: /s
        model.add_component(Parameter("kp10", 5.556e-06)) # binding of Sos to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km10", 0.06))      # dissociation of Sos from receptor, units: /s
        model.add_component(Parameter("kp11", 1.25e-06))  # binding of Grb2-Sos to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km11", 0.03))      # diss. of Grb2-Sos from receptor, units: /s
        model.add_component(Parameter("kp13", 2.5e-05))   # binding of Shc to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km13", 0.6))       # diss. of Shc from receptor, units: /s
        model.add_component(Parameter("kp15", 2.5e-07))   # binding of pShc to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km15", 0.3))       # diss. of pShc from receptor, units: /s
        model.add_component(Parameter("kp17", 1.667e-06)) # binding of Grb2 to RP-pShc (scaled), units: /molecule/s
        model.add_component(Parameter("km17", 0.1))       # diss. of Grb2 from RP-pShc, units: /s
        model.add_component(Parameter("kp18", 2.5e-07))   # binding of pShc-Grb2 to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km18", 0.3))       # diss. of pShc-Grb2 from receptor, units: /s
        model.add_component(Parameter("kp19", 5.556e-06)) # binding of Sos to RP-pShc-Grb2 (scaled), units: /molecule/s
        model.add_component(Parameter("km19", 0.0214))    # diss. of Sos from RP-pShc-Grb2, units: /s
        model.add_component(Parameter("kp20", 6.667e-08)) # binding of pShc-Grb2-Sos to receptor (scaled), units: /molecule/s
        model.add_component(Parameter("km20", 0.12))      # diss. of pShc-Grb2-Sos from receptor, units: /s
        model.add_component(Parameter("kp24", 5e-06))     # binding of Grb2-Sos to RP-pShc (scaled), units: /molecule/s
        model.add_component(Parameter("km24", 0.0429))    # diss. of Grb2-Sos from RP-pShc, units: /s
        model.add_component(Parameter("kp21", 1.667e-06)) # binding of pShc to Grb2 in cytosol (scaled), units: /molecule/s
        model.add_component(Parameter("km21", 0.01))      # diss. of Grb2 and SchP in cytosol, units: /s
        model.add_component(Parameter("kp23", 1.167e-05)) # binding of pShc to Grb2-Sos in cytosol (scaled), units: /molecule/s
        model.add_component(Parameter("km23", 0.1))       # diss. of Grb2-Sos and SchP in cytosol, units: /s
        model.add_component(Parameter("kp12", 5.556e-08)) # binding of Grb2 to Sos in cytosol (scaled), units: /molecule/s
        model.add_component(Parameter("km12", 0.0015))    # diss. of Grb2 and Sos in cytosol, units: /s
        model.add_component(Parameter("kp22", 1.667e-05)) # binding of pShc-Grb2 to Sos in cytosol (scaled), units: /molecule/s
        model.add_component(Parameter("km22", 0.064))     # diss. of pShc-Grb2 and Sos in cytosol, units: /s
        model.add_component(Parameter("kf_sos_ras", 1e-7)) # guess
        model.add_component(Parameter("kr_sos_ras", 1e-1)) # guess
        model.add_component(Parameter("kc_sos_ras", 1e-1)) # guess
        # Generic rates for MAPK cascade kinase/phosphatase binding, unbinding and catalysis.
        model.add_component(Parameter('kfc_raf_mek', 1e-6))
        model.add_component(Parameter('kfc_mek_erk', 1e-7))
        model.add_component(Parameter('kr_bind', 1e-1))
        model.add_component(Parameter('kcat_phos', 5.0))
        model.add_component(Parameter('kfc_dephos', 1e-6))

        model.rules['EGFR_EGF_bind'].rate_forward = model.parameters['kp1']
        model.rules['EGFR_EGF_bind'].rate_reverse = model.parameters['km1']
        model.rules['EGFR_EGFR_bind'].rate_forward = model.parameters['kp2']
        model.rules['EGFR_EGFR_bind'].rate_reverse = model.parameters['km2']
        model.rules['GRB2_EGFR_bind'].rate_forward = model.parameters['kp9']
        model.rules['GRB2_EGFR_bind'].rate_reverse = model.parameters['km9']
        model.rules['GRB2_SOS1_bind'].rate_forward = model.parameters['kp12']
        model.rules['GRB2_SOS1_bind'].rate_reverse = model.parameters['km12']
        model.rules['SOS1_HRAS_bind'].rate_forward = model.parameters['kf_sos_ras']
        model.rules['SOS1_HRAS_bind'].rate_reverse = model.parameters['kr_sos_ras']
        model.rules['HRAS_GTP_bind'].rate_forward = model.parameters['kf_bind']
        model.rules['HRAS_GTP_bind'].rate_reverse = model.parameters['kr_bind']
        model.rules['HRAS_RAF1_bind'].rate_forward = model.parameters['kf_bind']
        model.rules['HRAS_RAF1_bind'].rate_reverse = model.parameters['kr_bind']
        model.rules['RAF1_autophospho_RAF1_T491'].rate_forward = model.parameters['kf_autophospho']
        model.rules['RAF1_autophospho_RAF1_S494'].rate_forward = model.parameters['kf_autophospho']
        model.rules['RAF1_phospho_MAP2K1_S218'].rate_forward = model.parameters['kfc_raf_mek']
        model.rules['RAF1_phospho_MAP2K1_S222'].rate_forward = model.parameters['kfc_raf_mek']
        model.rules['MAP2K1_phospho_MAPK1_Y187'].rate_forward = model.parameters['kfc_mek_erk']
        model.rules['MAP2K1_phospho_MAPK1_T185'].rate_forward = model.parameters['kfc_mek_erk']
        model.rules['DUSP3_dephospho_MAPK1_T185'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP4_dephospho_MAPK1_T185'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP6_dephospho_MAPK1_T185'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP7_dephospho_MAPK1_T185'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP3_dephospho_MAPK1_Y187'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP4_dephospho_MAPK1_Y187'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP6_dephospho_MAPK1_Y187'].rate_forward = model.parameters['kfc_dephos']
        model.rules['DUSP7_dephospho_MAPK1_Y187'].rate_forward = model.parameters['kfc_dephos']


        # Set initial conditions
        set_initial(model, 'BT20')
        model.parameters['EGF_0'].value = 1e6

        # Set observables
        print 'Setting observables...'
        SOS1 = model.monomers['SOS1']
        HRAS = model.monomers['HRAS']
        RAF1 = model.monomers['RAF1']
        MAP2K1 = model.monomers['MAP2K1']
        MAPK1 = model.monomers['MAPK1']
        model.add_component(Observable('SOSbound', SOS1(GRB2=ANY)))
        model.add_component(Observable('RASGTP', HRAS(GTP=ANY)))
        model.add_component(Observable('RAFPP', RAF1(T491='p', S494='p')))
        model.add_component(Observable('MEKPP', MAP2K1(S218='p', S222='p')))
        model.add_component(Observable('ERKPP', MAPK1(T185='p', Y187='p')))

        # Set EGF dose
        model.parameters['EGF_0'].value = 1e6

        # Run model simulation
        # TODO: save figure as file
        ts = numpy.linspace(0, 120*60, 1000)

        import pickle
        if args.loadmodel:
            print 'Loading ODE solver...'
            with open('model_odes.pkl', 'rb') as fh:
                model = pickle.load(fh)
            solver = Solver(model, ts)
        else:
            print 'Constructing ODE solver...'
            tstart = time.time()
            solver = Solver(model, ts)
            tend = time.time()
            print '> Solver construction took %ds' % (tend - tstart)
            with open('model_odes.pkl', 'wb') as fh:
                pickle.dump(model, fh)

        print 'Running simulation...'
        tstart = time.time()
        solver.run()
        tend = time.time()
        print '> Simulation took %ds' % (tend - tstart)
        plt.ion()
        # plt.plot(ts, solver.yobs['SOSbound'] / model.parameters['SOS1_0'].value, label='SOS-bound')
        # plt.plot(ts, solver.yobs['RASGTP'] / model.parameters['HRAS_0'].value, label='HRAS-GTP')
        plt.plot(ts, solver.yobs['RAFPP'] / model.parameters['RAF1_0'].value, label='RAF1-pTpS')
        plt.plot(ts, solver.yobs['MEKPP'] / model.parameters['MAP2K1_0'].value, label='MEK1-pSpS')
        plt.plot(ts, solver.yobs['ERKPP'] / model.parameters['MAPK1_0'].value, label='ERK2-pTpY')
        plt.xlabel('Time (s)', fontsize=18)
        plt.ylabel('Normalized amount', fontsize=18)
        plt.xlim([0, ts[-1]])
        plt.legend(loc='lower left', fontsize=14)
        plt.show()

        print 'Running dose-response simulation...'
        egf_dose = numpy.logspace(3, 9, 50)
        erk_response = []
        tstart = time.time()
        for ed in egf_dose:
            model.parameters['EGF_0'].value = ed
            solver.run()
            erk_response.append(solver.yobs['ERKPP'][500])
        tend = time.time()
        print '> Dose response simulation took %ds' % (tend - tstart)
        plt.figure()
        plt.plot(egf_dose, erk_response, label='ERK-pTpY response')
        plt.xscale('log')
        plt.xlabel('EGF dose', fontsize=18)
        plt.ylabel('ERK2-pTpY amount at 1hr', fontsize=18)
        plt.legend(loc='lower right', fontsize=14)
