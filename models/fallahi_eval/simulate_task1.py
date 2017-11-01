import copy
import numpy
from pysb.bng import generate_equations
from pysb.simulator import ScipyOdeSimulator
from pysb import Observable, ReactionPattern, ComplexPattern

from util import pkldump, pklload
from process_data import *

# ### User defined parameters
time_lim_hrs = 10
time_steps = 100
# ##########################


def run_simulation_condition(sim, inhib_monomers):
    if inhib_monomers:
        print('Setting inhibition condition for %s' %
              ', '.join(inhib_monomers))
    old_params = {p.name: p.value for p in sim.model.parameters}
    sim_params = copy.copy(old_params)
    for monomer_name in inhib_monomers:
        init_name = monomer_name + '_0'
        sim_params[init_name] = 0

    print('Running simulation')
    res = sim.run(param_values=sim_params)

    print('Reverting change to model initial conditions')
    for monomer_name in inhib_monomers:
        init_name = monomer_name + '_0'
        model.parameters[init_name].value = old_params[init_name]
    return res


def set_model_observables(model):
    print('Setting model observables')
    s6 = model.monomers['RPS6']
    site_patterns = [{'S235': 'p', 'S236': 'u'},
                     {'S235': 'u', 'S236': 'p'},
                     {'S235': 'p', 'S236': 'p'}]
    pattern_terms = [ComplexPattern([s6(**pat)], None)
                     for pat in site_patterns]
    obs_pattern = ReactionPattern(pattern_terms)
    obs = Observable('pS6', obs_pattern, _export=False)
    model.add_component(obs)


if __name__ == '__main__':
    time_lim_seconds = time_lim_hrs * 3600
    ts = numpy.linspace(0, time_lim_seconds, time_steps + 1)

    for cell_line  in ('C32', 'LOXIMVI', 'MMACSF', 'MZ7MEL', 'RVH421'):
        print('Loading model for %s' % cell_line)
        model = pklload('pysb_model_%s' % cell_line)

        set_model_observables(model)

        print('Generating model equations')
        generate_equations(model)

        print('Instantiating simulator')
        sim = ScipyOdeSimulator(model, ts, use_theano=True, verbose=True)

        res_drug = run_simulation_condition(sim, ['ARAF', 'BRAF', 'RAF1'])
        res_nodrug = run_simulation_condition(sim, [])
        break

