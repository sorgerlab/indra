import numpy as np
import matplotlib.pyplot as plt
from pysb.integrate import Solver
import assemble_p53model


def run_model(model_name):
    model = assemble_p53model.assemble_model(model_name)
    sim_hours = 200
    ts = np.linspace(0, sim_hours*3600, sim_hours*60)
    solver = Solver(model, ts)
    solver.run()
    plt.plot(ts, solver.yobs['P53_active'], 'r')
    plt.xticks([])
    plt.xlabel('time (a.u)', fontsize=15)
    plt.ylabel('Active p53')
    plt.yticks([])
    plt.savefig(model_name)      
    return ts, solver


