import numpy as np
from pysb.integrate import Solver
from matplotlib import pyplot as plt
from indra import bel, biopax, trips
from indra.util import plot_formatting as pf
from indra.assemblers import PysbAssembler

# 1. TEXT
# User defines text:
text = ('MEK1 phosphorylates ERK2 on threonine 185 and tyrosine 187.')

# Show round trip going out to TRIPS/DRUM web service,
# return logical form to INDRA, which is queried by
# INDRA for relevant statements
tp = trips.process_text(text)

# Now generate PySB model
stmts = tp.statements # Don't show this one

pa = PysbAssembler()
pa.add_statements(stmts)
pa.make_model()
t = np.linspace(0, 2500)
sol = Solver(pa.model, t)
sol.run()

pf.set_fig_params()
plt.ion()
plt.figure(figsize=(1, 1), dpi=300)
species_names = [str(s) for s in pa.model.species]
plt.plot(t, sol.y[:, species_names.index("MAPK1(T185='u', Y187='u')")],
         'b', label='MAPK1.uu')
plt.plot(t, sol.y[:, species_names.index("MAPK1(T185='p', Y187='u')")] +
            sol.y[:, species_names.index("MAPK1(T185='u', Y187='p')")],
         'g', label='MAPK1.p')
plt.plot(t, sol.y[:, species_names.index("MAPK1(T185='p', Y187='p')")],
         'r', label='MAPK1.pp')
plt.xlabel('Time')
plt.ylabel('Amount')
plt.legend(loc='upper right', fontsize=4)
plt.xticks([])
plt.yticks([])
pf.format_axis(plt.gca())
