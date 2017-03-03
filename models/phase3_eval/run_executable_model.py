from copy import copy, deepcopy
import numpy
import pickle
import matplotlib.pyplot as plt
from pysb.integrate import ScipyOdeSimulator

from process_data import get_drug_targets

ts = numpy.linspace(0, 86400, 1440)

def load_model(pkl_file):
    with open(pkl_file, 'r') as fh:
        model = pickle.load(fh)
        return model

def sim_inhibit(sim, model, inhibit_names, inhibit_ratio):
    params = [p.value for p in model.parameters]
    for i, p in enumerate(model.parameters):
        for name in inhibit_names:
            pname = '%s_0' % name
            if p.name == pname:
                params[i] = params[i] * (1.0 - inhibit_ratio)

    print('Running perturbed condition')
    sim._initials = None
    res = sim.run(param_values=params)
    res = deepcopy(res.observables)
    return res

def sim_unperturbed(sim):
    print('Running unperturbed condition')
    sim._initials = None
    params = [p.value for p in model.parameters]
    res = sim.run(param_values=params)
    res = deepcopy(res.observables)
    return res

def plot_all(model, all_results):
    for obs in model.observables:
        fname = 'output/%s.png' % obs.name
        plt.figure()
        for drug_abbrev in all_results.keys():
            plt.plot(ts[:60], all_results[drug_abbrev][obs.name][:60],
                     label=drug_abbrev)
        plt.xlabel('Time (seconds)')
        plt.ylabel('Amount (molecules)')
        plt.legend()
        plt.savefig(fname)

def compare_all(model, all_results):
    ratios = {}
    for obs in model.observables:
        ratios[obs.name] = {}
        base = all_results['unperturbed'][obs.name][10]
        for drug_abbrev in all_results.keys():
            if drug_abbrev == 'unperturbed':
                continue
            treated = all_results[drug_abbrev][obs.name][10]
            if base < 1e-10:
                ratio = 1.0
            else:
                ratio = treated / base
            ratios[obs.name][drug_abbrev] = ratio
    return ratios

def sim_all(model, drug_targets):
    all_results = {}
    sim = ScipyOdeSimulator(model, ts)
    res = sim_unperturbed(sim)
    all_results['unperturbed'] = res
    for drug_abbrev, targets in drug_targets.items():
        res = sim_inhibit(sim, model, targets, 0.99)
        all_results[drug_abbrev] = res
    with open('output/sim_results.pkl', 'wb') as fh:
        pickle.dump(all_results, fh)
    return all_results

if __name__ == '__main__':
    print('Loading model')
    model = load_model('output/korkut_model_pysb_odes.pkl')
    drug_targets = get_drug_targets('data/drug_grounding.csv')
    rerun = True
    if rerun:
        all_results = sim_all(model, drug_targets)
    else:
        with open('output/sim_results.pkl', 'r') as fh:
            all_results = pickle.load(fh)
    plot_all(model, all_results)
