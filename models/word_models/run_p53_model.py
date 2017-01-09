import os
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
            tp = trips.process_text(fh.read(), xml_fname)

    pa = PysbAssembler()
    pa.add_statements(tp.statements)
    model = pa.make_model()

    p53 = model.monomers['TP53']
    obs = Observable('P53_active', p53(activity='active'))
    model.add_component(obs)
    if not model_name.endswith('var'):
        model.parameters['kf_aa_act_1'].value = 5e-06
    model.parameters['kf_pt_act_1'].value = 1e-05

    if model_name == 'p53_ATM':
        model.add_component(Parameter('ATMa_0', 1))
        atm = model.monomers['ATM']
        model.initial(atm(activity='active'),
                      model.parameters['ATMa_0'])
        model.parameters['kf_pa_act_1'].value = 1e-04

    if model_name == 'p53_ATR':
        model.add_component(Parameter('ATRa_0', 1))
        atr = model.monomers['ATR']
        model.initial(atr(activity='active'),
                      model.parameters['ATRa_0'])

    pa.model = model
    pa.save_model('%s.py' % model_name)
    return model

def run_model(model):
    sim_hours = 200
    ts = np.linspace(0, sim_hours*3600, sim_hours*60)
    solver = Solver(model, ts)
    solver.run()
    plt.figure(figsize=(2,2), dpi=300)
    set_fig_params()
    plt.plot(ts, solver.yobs['P53_active'], 'r')
    plt.xticks([])
    plt.xlabel('Time (a.u.)', fontsize=15)
    plt.ylabel('Active p53')
    plt.yticks([])
    plt.savefig(model_name + '.pdf')
    return ts, solver

if __name__ == '__main__':
    model_names = ['p53_ATR', 'p53_ATM', 'p53_ATM_var']
    #model_names = ['p53_ATM_var']
    for model_name in model_names:
        model = assemble_model(model_name, reread=False)
        ts, solver = run_model(model)
