import numpy as np
from pysb.integrate import Solver
from pysb.export import export
from matplotlib import pyplot as plt
from indra.sources import biopax
from indra.sources import bel, trips
from indra.util import plot_formatting as pf
from indra.assemblers import PysbAssembler

def plot_result(model, sol):
    pf.set_fig_params()
    plt.ion()
    plt.figure(figsize=(1, 1), dpi=300)
    species_names = [str(s) for s in model.species]
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


def export_hello(model, formats):
    for f in formats:
        model_export = export(model, f)
        extension = (f if f != 'pysb_flat' else 'py')
        fname = 'hello_indra_model.%s' % extension
        with open(fname, 'wb') as fh:
            fh.write(model_export)


# User defines text
text = 'MEK1 phosphorylates ERK2 on threonine 185 and tyrosine 187.'

# Process text using TRIPS processor
tp = trips.process_text(text)

# Get the list of extracted Statements
stmts = tp.statements

# Assemble a PySB model
pa = PysbAssembler()
pa.add_statements(stmts)
pa.make_model()

# Run simulation
t = np.linspace(0, 300)
sol = Solver(pa.model, t)
sol.run()

# Plot the result
plot_result(pa.model, sol)

# Export model
export_hello(pa.model, ['sbml', 'bngl', 'kappa', 'pysb_flat'])
