from indra import bel, biopax, trips
from indra.assemblers import PysbAssembler
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt

# 1. TEXT
# User defines text:
text = ('MEK1 phosphorylates ERK2 on Thr-185 and Tyr-187. DUSP4 '
       'dephosphorylates ERK2 on Thr-185 and Tyr-187.')

# Show round trip going out to TRIPS/DRUM web service,
# return logical form to INDRA, which is queried by
# INDRA for relevant statements
tp = trips.process_text(text)

# 2. BIOPAX
# Show round trip going out to PC web service, returning
# Biopax model, queried by INDRA
#bp = biopax.process_pc_pathsbetween(['MAP2K1', 'MAPK1'])
#bp.get_phosphorylation()

# 3. BEL
# Show round trip to NDeX/BEL service
#belp = bel.process_ndex_neighborhood(['MAP2K1'])

# Statements can now be collected since they are in the same format
#stmts = tp.statements + bp.statements + belp.statements

# .... (do stuff with statements here, filter, curate, etc.)

# Now generate PySB model
stmts = tp.statements # Don't show this one

pa = PysbAssembler()
pa.add_statements(stmts)
pa.make_model()
t = np.linspace(0, 10000)
sol = Solver(pa.model, t)
sol.run()

labels = [str(s) for s in pa.model.species]
plt.ion()
for i in range(sol.y.shape[1]):
    plt.plot(t, sol.y[:, i], label=labels[i])
plt.legend(loc='upper right')
