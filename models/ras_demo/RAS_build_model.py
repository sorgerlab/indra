import pickle
import numpy as np
import sys
from pysb import Parameter, Observable, Rule
from pysb.integrate import Solver
import pysb.core

from indra.assemblers import PysbAssembler
from indra import bel

def get_base_state(mon):
    sdict = {}
    for s in mon.sites:
        if s in mon.site_states:
            sdict[s]=mon.site_states[s][0]
        else:
            sdict[s]=None
    return mon(sdict)

def add_initial(model, pattern, value):
    """Add initial condition; if an initial condition for this pattern
    already exists, override it."""
    complex_pattern = pysb.as_complex_pattern(pattern)

    for other_cp, other_value in model.initial_conditions:
        if complex_pattern.is_equivalent_to(other_cp):
            model.initial_conditions.remove((other_cp, other_value))
    model.initial(complex_pattern, value)

# Generate model rules via indra
pa = PysbAssembler()
bp = bel.process_belrdf('../../data/RAS_combined.rdf')
pa.add_statements(bp.statements)
model = pa.make_model(initial_conditions=False)

# Useful shortcuts to access model components
m = model.monomers
p = model.parameters

# Add ligand to the model
model.add_component(Monomer('EGF'))
model.add_component(Parameter('kf_ee_act', 1e-6))
model.add_component(
    Rule('EGF_activates_EGFR',
         m['EGF']() + m['EGFR']({'Kinase':'inactive'}) >>
         m['EGF']() + m['EGFR']({'Kinase':'active'}),
         p['kf_ee_act']))

# Add initial conditions
# Values taken from Chen et al. (2009) Mol. Syst. Biol.
model.initial_conditions = []
model.add_component(Parameter('init_default',10000))
model.add_component(Parameter('EGF_0',10000))
model.add_component(Parameter('EGFR_0',70000))
model.add_component(Parameter('SOS_0',5000000))
model.add_component(Parameter('RAS_0',58000))
model.add_component(Parameter('RAF_0',40000))
model.add_component(Parameter('MAP2K1_0',3020000))
model.add_component(Parameter('MAPK1_0',695000))
model.add_component(Parameter('AKT_0',905000))

add_initial(model, m['EGF'](), p['EGF_0'])
add_initial(model, get_base_state(m['EGFR']), p['EGFR_0'])
add_initial(model, get_base_state(m['SOS1']), p['SOS_0'])
add_initial(model, get_base_state(m['HRAS']), p['RAS_0'])
add_initial(model, get_base_state(m['NRAS']), p['RAS_0'])
add_initial(model, get_base_state(m['KRAS']), p['RAS_0'])
add_initial(model, get_base_state(m['ARAF']), p['RAF_0'])
add_initial(model, get_base_state(m['BRAF']), p['RAF_0'])
add_initial(model, get_base_state(m['RAF1']), p['RAF_0'])
add_initial(model, get_base_state(m['MAP2K1']), p['MAP2K1_0'])
add_initial(model, get_base_state(m['MAPK1']), p['MAPK1_0'])
add_initial(model, m['RASA2'](Catalytic='active'), p['init_default'])
add_initial(model, m['RASA3'](Catalytic='active'), p['init_default'])

# Add observables
model.add_component(Observable("ERKact", m['MAPK1'](Kinase='active')))
model.add_component(Observable("MEKact", m['MAP2K1'](Kinase='active')))

# Import solver and generate model equations
t = np.linspace(0,25,11)
solver = Solver(model,t)

# Pickle the model
with open('RAS_combined_model.pkl','wb') as fh:
    pickle.dump(model,fh)
