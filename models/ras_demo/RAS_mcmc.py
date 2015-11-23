from pysb.integrate import Solver
import pickle
import numpy as np
import scipy
import matplotlib.pyplot as plt
import bayessb
import time
import sys


def likelihood(mcmc,position):
    yobs = mcmc.simulate(position,observables=True)
    ll = np.sum((yobs['ERKact'][1:]-mcmc.options.exp_data) ** 2 /
                (2.0 * mcmc.options.exp_data_std**2))
    return ll

def prior(mcmc,position):
    lp = np.sum((position - mcmc.options.prior_mean) ** 2 /
                (2.0 * mcmc.options.prior_var))
    return lp

def step(mcmc):
    """Print out some statistics every 20 steps"""
    if mcmc.iter % 20 == 0:
        print 'iter=%-5d  sigma=%-.3f  T=%-.3f  acc=%-.3f, lkl=%g  ' \
              'prior=%g  post=%g' % \
              (mcmc.iter, mcmc.sig_value, mcmc.T,
               float(mcmc.acceptance)/(mcmc.iter+1), mcmc.accept_likelihood,
               mcmc.accept_prior, mcmc.accept_posterior)

class Object:
    def __init__(self):
        self.model              = None
        self.estimate_params    = None
        self.initial_values     = None
        self.tspan              = None
        self.step_fn            = None
        self.likelihood_fn      = None
        self.prior_fn           = None
        self.nsteps             = None
        self.use_hessian        = False
        self.start_random       = False
        self.boundary_option    = False
        self.rtol               = None
        self.atol               = None
        self.norm_step_size     = 0.75
        self.hessian_period     = 25000
        self.hessian_scale      = 0.085
        self.sigma_adj_interval = None
        self.anneal_length      = None
        self.T_init             = 10
        self.accept_rate_target = 0.3
        self.sigma_max          = 1
        self.sigma_min          = 0.25
        self.sigma_step         = 0.125
        self.thermo_temp        = 1
        self.seed               = None

    def copy(self):
        new_options = Object()
        new_options.__dict__.update(self.__dict__)
        return new_options

model_fname = 'RAS_combined_model.pkl'
try:
    model = pickle.load(open(model_fname,'rb'))
except IOError:
    print 'Could not open model file %s' % model_fname
    sys.exit()

model.integrator = scipy.integrate.ode(model.odes)
model.integrator.set_integrator('vode', method='bdf', with_jacobian=True,
                                rtol=1e-3,atol=1e-6,nsteps=20000,order=5)

opts = Object()

tplot = np.linspace(0,25,101)
tdata = np.linspace(0,25,11)

solver = Solver(model,tdata)
solver.run()
sim_data = solver.yobs['ERKact'][1:]

opts.exp_data_std = 40000
opts.exp_data = sim_data + np.random.randn(len(sim_data))*opts.exp_data_std
opts.nsteps = 10000
opts.model = model
opts.tspan = tdata
opts.likelihood_fn = likelihood
opts.prior_fn = prior
opts.step_fn = step
opts.estimate_params = map(model.parameters.get,
                           ['kf_ee_act', 'kf_bh_bind_1', 'kf_sn_act_1',
                            'kf_ha_act_1'])
opts.prior_mean = [np.log10(p.value) for p in opts.estimate_params]
opts.prior_var = np.empty_like(opts.prior_mean)
opts.prior_var.fill(1.0)
opts.seed = 123
opts.use_hessian = True
opts.anneal_length = 200
opts.norm_step_size = 0.05

mcmc = bayessb.MCMC(opts)

start_time = time.time()
mcmc.run()
end_time = time.time()

print "Elapsed time", end_time - start_time, "s"

mcmc.options.tspan = tplot
mcmc.solver.__init__(mcmc.solver.model, mcmc.options.tspan)

plt.ion()
for p in mcmc.positions[::10]:
    yobs = mcmc.simulate(p,observables=True)
    plt.plot(mcmc.options.tspan,yobs['ERKact'],'b-',alpha=0.1)

plt.plot(0, 0, 'b-', label='Simulation')
plt.errorbar(np.linspace(0,25,11)[1:], opts.exp_data, opts.exp_data_std,
             fmt='ro',label='Data')
plt.ylim(ymin=-50000)
plt.legend(loc='upper left')
plt.show()
