from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from indra import trips
from indra.assemblers import PysbAssembler
from indra.util.plot_formatting import *
from pysb import Observable, Parameter
from pysb.integrate import Solver

def assemble_model(model_name, reread=False):
    xml_fname = model_name + '.xml'
    if not reread:
        print('Processing %s' % xml_fname)
        if os.path.exists(xml_fname):
            with open(xml_fname, 'rb') as fh:
                tp = trips.process_xml(fh.read())
        else:
            reread = True
    if reread:
        fname = model_name + '.txt'
        print('Reading %s' % fname)
        with open(fname, 'rb') as fh:
            ts = time.time()
            tp = trips.process_text(fh.read(), xml_fname)
            te = time.time()
            print('Reading took %.2fs' % (te-ts))
    print('Assembling statements:')
    for i, st in enumerate(tp.statements):
        print('%d: %s' % (i, st))
    print('----------------------')

    pa = PysbAssembler()
    pa.add_statements(tp.statements)
    ts = time.time()
    model = pa.make_model()
    te = time.time()
    print('Assembly took %.2fs' % (te-ts))
    model.name = model_name

    p53 = model.monomers['TP53']
    obs = Observable(b'p53_active', p53(activity='active'))
    model.add_component(obs)
    if not model_name.endswith('var'):
        model.parameters['kf_aa_act_1'].value = 5e-06
    model.parameters['kf_pt_act_1'].value = 5e-06

    if model_name == 'p53_ATM':
        model.add_component(Parameter('ATMa_0', 1))
        atm = model.monomers['ATM']
        model.initial(atm(activity='active'),
                      model.parameters['ATMa_0'])
        model.parameters['kf_pa_act_1'].value = 1e-04
        obs = Observable(b'atm_active', atm(activity='active'))
        model.add_component(obs)

    if model_name == 'p53_ATR':
        model.add_component(Parameter('ATRa_0', 1))
        atr = model.monomers['ATR']
        model.initial(atr(activity='active'),
                      model.parameters['ATRa_0'])
        obs = Observable(b'atr_active', atr(activity='active'))
        model.add_component(obs)

    if model_name == 'p53_ATM_var':
        model.add_component(Parameter('ATMa_0', 1))
        atm = model.monomers['ATM']
        model.initial(atm(phospho='p'),
                      model.parameters['ATMa_0'])
        model.parameters['kf_pa_dephosphorylation_1'].value = 1e-04
        model.parameters['MDM2_0'].value = 0
        model.parameters['kf_m_deg_1'].value = 8e-01
        model.parameters['kf_tm_synth_1'].value = 0.2
        model.parameters['kf_aa_phosphorylation_1'].value = 5e-06
        obs = Observable(b'atm_active', atm(phospho='p'))
        model.add_component(obs)

    pa.model = model
    pa.save_model('%s.py' % model_name)
    return model

def run_model(model):
    sim_hours = 20
    ts = np.linspace(0, sim_hours*3600, sim_hours*60)
    tst = time.time()
    solver = Solver(model, ts)
    solver.run()
    te = time.time()
    print('Simulation took %.2fs' % (te-tst))
    plt.figure(figsize=(2,2), dpi=300)
    set_fig_params()
    plt.plot(ts, solver.yobs['p53_active'], 'r')
    plt.xticks([])
    plt.xlabel('Time (a.u.)', fontsize=12)
    plt.ylabel('Active p53', fontsize=12)
    plt.yticks([])
    plt.savefig(model.name + '.pdf')
    return ts, solver

if __name__ == '__main__':
    reread = False
    model_names = ['p53_ATR', 'p53_ATM', 'p53_ATM_var']
    for model_name in model_names:
        model = assemble_model(model_name, reread=reread)
        ts, solver = run_model(model)
