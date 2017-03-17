import numpy
import matplotlib.pyplot as plt
from indra.util import plot_formatting as pf

# Import ATM models
from ATM_v1 import model as ATM_v1
from ATM_v2 import model as ATM_v2
from ATM_v3 import model as ATM_v3
from ATM_v4a import model as ATM_v4a
from ATM_v4b import model as ATM_v4b

# Import ATR models
from ATR_v1 import model as ATR_v1
from ATR_v2 import model as ATR_v2
from ATR_v3 import model as ATR_v3

from run_p53_model import run_model
from pysb.bng import generate_equations
from pysb.simulator import ScipyOdeSimulator as Solver

def sample_params(mu, sigma):
    r = numpy.random.randn(len(mu))
    p = numpy.power(10, mu + r*sigma)
    return p

def parameter_sweep(model, sigma, ns):
    generate_equations(model)
    logp = [numpy.log10(p.value) for p in model.parameters]
    ts = numpy.linspace(0, 20*3600, 20*60)
    solver = Solver(model, ts)
    pf.set_fig_params()
    plt.figure(figsize=(1.8, 1), dpi=300)
    for i in range(ns):
        psample = sample_params(logp, 0.05)
        res = solver.run(param_values=psample)
        signal = res.observables['p53_active']
        plt.plot(signal, color=(0.7, 0.7, 0.7), alpha=0.3)
    # Highlighted
    colors = ['g', 'y', 'c']
    for c in colors:
        psample = sample_params(logp, 0.05)
        res = solver.run(param_values=psample)
        signal = res.observables['p53_active']
        plt.plot(signal, c)

    # Nominal
    solver = Solver(model, ts)
    res = solver.run()
    signal = res.observables['p53_active']
    plt.plot(signal, 'r')

    plt.xticks([])
    plt.xlabel('Time (a.u.)', fontsize=7)
    plt.ylabel('Active p53', fontsize=7)
    plt.yticks([])
    plt.ylim(ymin=0)
    pf.format_axis(plt.gca())
    plt.savefig(model.name + '_sample.pdf')

if __name__ == '__main__':
    models = [ATM_v1, ATM_v2, ATM_v3, ATM_v4a, ATM_v4b,
              ATR_v1, ATR_v2, ATR_v3]
    for model in models:
        print(model.name)
        parameter_sweep(model, 1, 50)

