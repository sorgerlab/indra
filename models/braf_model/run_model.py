from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import pickle
import matplotlib
import numpy as np
# To suppress plots from popping up
import matplotlib.pyplot as plt
from pysb.integrate import Solver
from pysb.bng import generate_equations
from pysb.export import export
from indra.util.plot_formatting import *

import assemble_model


def simulate_untreated(model, ts, egf_dose=1000):
    model.parameters['EGF_0'].value = egf_dose
    model.parameters['VEMURAFENIB_0'].value = 0
    solver = Solver(model, ts)
    solver.run()
    return (solver.yobs, solver.y)


def simulate_vemurafenib_treatment(model, ts, y0):
    # New initial conditions for simulation post event
    y_init = y0
    y_init[model.observables['Vem_free'].species[0]] = 2e5

    # Continue model simulation with y_init as new initial condition
    solver = Solver(model, ts)
    solver.run(y0=y_init)

    # Concatenate the two simulations
    return solver.yobs, solver.y


def is_steady_state(y):
    for k in y.dtype.fields.keys():
        d = np.diff(y[k])
        if not all(abs(d) < max(abs(y[k]))*0.01):
            return False
    return True


def get_steady_state(model, y0):
    sim_hours = 1
    t = np.linspace(0, sim_hours*3600, sim_hours*60)
    solver = Solver(model, t)
    ss = False
    y_init = y0
    while not ss:
        solver.run(y0=y_init)
        ss = is_steady_state(solver.yobs)
        y_init = solver.y[-1]
    steady_state = solver.yobs[-1]
    timecourse = solver.yobs
    return steady_state, timecourse


def plot_fold_change_time(ts, yobs, yobs_ref, save_plot):
    erk_foldchange = yobs['ERK_p'] / yobs_ref['ERK_p']
    ras_foldchange = yobs['RAS_active'] / yobs_ref['RAS_active']
    plt.figure(figsize=(2, 2), dpi=300)
    set_fig_params()
    plt.plot(ts, erk_foldchange, 'b')
    plt.plot(ts, ras_foldchange, 'g')
    plt.xticks([])
    plt.xlabel('Time (a.u)', fontsize=7)

    plt.xlim([0, 30000])
    plt.ylabel('Fold-change after Vemurafenib treatment', fontsize=7)
    plt.legend(['Phosphorylated ERK', 'Active RAS'], fontsize=7)
    ax = plt.gca()
    format_axis(ax)
    plt.savefig(save_plot)
    plt.clf()


def get_egf_vem_doseresp(egf_doses, vem_doses, readout='absolute'):
    erk = np.zeros((len(egf_doses), len(vem_doses)))
    ras = np.zeros((len(egf_doses), len(vem_doses)))
    for i, ed in enumerate(egf_doses):
        yobs0, y0 = simulate_untreated(model, ts, ed)
        for j, vd in enumerate(vem_doses):
            y0_last = y0[-1]
            y0_last[model.observables['Vem_free'].species[0]] = vd
            yobs_ss, yobs_tc = get_steady_state(model, y0_last)
            if readout == 'absolute':
                erk[i,j] = yobs_ss['ERK_p']
                ras[i,j] = yobs_ss['RAS_active']
            elif readout == 'foldchange':
                erk[i,j] = yobs_ss['ERK_p'] / yobs0['ERK_p'][-1]
                ras[i,j] = yobs_ss['RAS_active'] / yobs0['RAS_active'][-1]
            elif readout == 'ss_min_diff':
                erk[i,j] = yobs_ss['ERK_p'] - yobs_tc['ERK_p'].min()
                ras[i,j] = yobs_ss['RAS_active']
    return erk, ras


def plot_egf_vem_dose(egf_doses, vem_doses, erk, ras, save_plot):
    cmap = plt.get_cmap('jet')
    plt.figure(figsize=(3, 3), dpi=300)
    set_fig_params()
    f, (ax1, ax2) = plt.subplots(1, 2)
    max_erk = np.amax(erk) + 1
    max_ras = np.amax(ras) + 1

    f1 = ax1.imshow(np.flipud(ras), cmap=cmap, vmin=0, vmax=max_ras,
                    interpolation='hermite')
    ax1.set_title('Active RAS')
    ax1.set_xlabel('[Vemurafenib]', fontsize=7)
    ax1.set_ylabel('[EGF]', fontsize=7)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_aspect('equal')
    cb = f.colorbar(f1, ax=ax1, fraction=0.05, pad=0.05)
    cb.set_ticks([0, max_ras])

    f2 = ax2.imshow(np.flipud(erk), cmap=cmap, vmin=0, vmax=max_erk,
                    interpolation='hermite')
    ax2.set_title('pERK rebound')
    ax2.set_xlabel('[Vemurafenib]', fontsize=7)
    ax2.set_ylabel('[EGF]', fontsize=7)
    ax2.set_aspect('equal')
    ax2.set_xticks([])
    ax2.set_yticks([])
    cb = f.colorbar(f2, ax=ax2, fraction=0.05, pad=0.05)
    cb.set_ticks([0, max_erk])

    format_axis(ax1)
    format_axis(ax2)
    f.savefig(save_plot)


def load_model(model_id):
    model = pickle.load(open('model%d.pkl' % model_id, 'rb'))
    return model


def save_model(model):
    with open('model%d.pkl' % model_id, 'wb') as fh:
        pickle.dump(model, fh)


def print_statistics(model):
    print('# monomers: %d' % len(model.monomers))
    print('# rules: %d' % len(model.rules))
    p_rate = [p for p in model.parameters if p.name.startswith('k')]
    print('# rate constants: %d' % len(p_rate))
    print('# initial conditions: %d' % len(model.initial_conditions))
    print('# ODEs: %d' % len(model.odes))


def export_memi(model, formats, version):
    for f in formats:
        model_export = export(model, f)
        extension = (f if f != 'pysb_flat' else 'py')
        fname = 'MEMI%s.%s' % (version, extension)
        with open(fname, 'wb') as fh:
            fh.write(model_export)


if __name__ == '__main__':
    sim_hours = 10
    ts = np.linspace(0, sim_hours*3600, sim_hours*60)
    egf_doses = np.logspace(1, 4, 9)
    vem_doses = np.logspace(4, 6, 9)
    for model_id in (1, 2, 3):
        print('Running model %d' % model_id)
        print('----------------')
        if os.path.exists('model%d.pkl' % model_id):
            model = load_model(model_id)
        else:
            model = assemble_model.assemble_model(model_id)
            generate_equations(model)
            save_model(model)
        print_statistics(model)

        # Run time-course simulation
        yobs1, y1 = simulate_untreated(model, ts)
        yobs2, y2 = simulate_vemurafenib_treatment(model, ts, y1[-1])
        plot_fold_change_time(ts, yobs2, yobs1[-1],
                              'model%s_vem_time.pdf' % model_id)

        # Run dose response experiments
        erk_foldchange, ras_foldchange = \
            get_egf_vem_doseresp(egf_doses, vem_doses, 'ss_min_diff')
        plot_egf_vem_dose(egf_doses, vem_doses,
                              erk_foldchange, ras_foldchange,
                              'model%s_egf_vem_dose.pdf' % model_id)
        # Save results of dose response runs
        with open('doseresp_result_model%d.pkl' % model_id, 'wb') as fh:
            pickle.dump([egf_doses, vem_doses, erk_foldchange, ras_foldchange], 
                        fh)

        # Export models in various formats
        version = '1.%d' % (model_id-1)
        formats = ['sbml', 'bngl', 'kappa', 'pysb_flat']
        export_memi(model, formats, version)
