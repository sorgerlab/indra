import numpy
from pysb.bng import generate_equations
from pysb.simulator import ScipyOdeSimulator
from pysb import Observable, ReactionPattern

from util import pkldump, pklload
from process_data import *

# ### User defined parameters
time_lim_hrs = 10
time_steps = 100
# ##########################


def run_simulation_condition(model, ts, inhib_monomers):
    if inhib_monomers:
        print('Setting inhibition condition for %s' %
              ', '.join(inhib_monomers))
    old_init = {}
    for monomer_name in inhib_monomers:
        init_name = monome_name + '_0'
        old_init[monomer_name] = model.parameters[init_name].value
        model.parameters[init_name].value = 0

    print('Instantiating simulator')
    sim = ScipyOdeSimulator(model, ts)

    print('Running simulation')
    res = sim.run()

    print('Reverting change to model initial conditions')
    for monomer_name in inhib_monomers:
        init_name = monome_name + '_0'
        model.parameters[init_name].value = old_init[monomer_name]
    return res


def set_model_observables(model):
    s6 = model.monomers['RPS6']
    site_patterns = [
        {'S235': 'p', 'S236': 'u'}
        {'S235': 'u', 'S236': 'p'}
        {'S235': 'p', 'S236': 'p'}
        ]
    pattern_terms = [s6(**pat) for pat in site_patterns]
    obs_pattern = ReactionPattern(pattern_terms)
    obs = Observable('pS6', obs_pattern)
    model.add_component(obs)


if __name__ == '__main__':
    time_lim_seconds = time_lim_hrs * 3600
    ts = numpy.linspace(0, time_lim_seconds, time_steps + 1)

    for cell_line  in ('C32', 'LOXIMVI', 'MMACSF', 'MZ7MEL', 'RVH421'):
        print('Loading model for %s' % cell_line)
        model = pklload('pysb_model_%s' % cell_line)
        print('Generating model equations for %s' % cell_line)
        generate_equations(model)
        res_nodrug = run_simulation_condition(model, ts, [])
        res_drug = run_simulation_condition(model, ts, ['ARAF', 'BRAF', 'RAF1'])
        break

