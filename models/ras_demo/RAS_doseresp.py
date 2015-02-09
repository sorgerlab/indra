from pysb.integrate import Solver
import pickle
import numpy as np
import scipy
import matplotlib.pyplot as plt
import time

#model = pickle.load(open('RAS_combined_model.pkl','rb'))
model.integrator = scipy.integrate.ode(model.odes)
model.integrator.set_integrator('vode', method='bdf', with_jacobian=True, rtol=1e-3,atol=1e-6,nsteps=20000,order=5)
 
t = np.linspace(0,30,101)
solver = Solver(model,t)

y = []
for EGF_0 in np.logspace(-2,6,9):
    model.parameters['EGF_0'].value = EGF_0
    solver.run()
    y.append(solver.yobs.copy())

plt.ion()
for yy in y:
	plt.plot(t,yy["ERKact"],'b-')
plt.title("active ERK")
