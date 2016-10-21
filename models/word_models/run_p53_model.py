import numpy as np
import matplotlib.pyplot as plt
from indra import trips
from indra.assemblers import PysbAssembler
from pysb import Observable, Parameter
from pysb.integrate import Solver

def assemble_model(model_name):
    tp = trips.process_text(open(model_name + '.txt').read())
    pa = PysbAssembler()
    pa.add_statements(tp.statements)
    model = pa.make_model()

    p53 = model.monomers['TP53']
    obs = Observable('P53_active', p53(act='active'))
    model.add_component(obs)
    model.parameters['kf_aa_act_1'].value = 5e-06
    model.parameters['kf_pt_act_1'].value = 1e-05
    model.parameters['CDKN2A_0'].value = 0
    model.parameters['PROTEASE_0'].value = 0
    protease = model.monomers['PROTEASE']
    model.add_component(Parameter('CDKN2A_act_0', 100))
    model.add_component(Parameter('PROTEASE_act_0', 100))
    model.initial(protease(act='active'),
                  model.parameters['PROTEASE_act_0'])
    cdkn2a = model.monomers['CDKN2A']
    model.initial(cdkn2a(act='active'),
                  model.parameters['CDKN2A_act_0'])

    if model_name == 'p53_ATM':
        model.add_component(Parameter('ATMa_0', 1))
        atm = model.monomers['ATM']
        model.initial(atm(act='active'),
                      model.parameters['ATMa_0'])
        model.parameters['kf_pa_act_1'].value = 1e-04
        model.parameters['ATM_0'].value = 99.0

    if model_name == 'p53_ATR':
        model.add_component(Parameter('ATRa_0', 1))
        atr = model.monomers['ATR']
        model.initial(atr(act='active'),
                      model.parameters['ATRa_0'])
        model.parameters['ATR_0'].value = 99.0

    pa.model = model
    pa.save_model('%s.py' % model_name)
    return model

def run_model(model):
    sim_hours = 200
    ts = np.linspace(0, sim_hours*3600, sim_hours*60)
    solver = Solver(model, ts)
    solver.run()
    plt.plot(ts, solver.yobs['P53_active'], 'r')
    plt.xticks([])
    plt.xlabel('time (a.u)', fontsize=15)
    plt.ylabel('Active p53')
    plt.yticks([])
    plt.savefig(model_name + '.png')
    return ts, solver

if __name__ == '__main__':
    model_name = 'p53_ATR'
    model = assemble_model(model_name)
    ts, solver = run_model(model)
