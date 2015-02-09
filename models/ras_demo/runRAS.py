import numpy as np
import sys
from pysb import *
sys.path.append('../../')
from belpy import rdf_to_pysb

def get_base_state(mon):
    sdict = {}
    for s in mon.sites:
        if mon.site_states.has_key(s):
            sdict[s]=mon.site_states[s][0]
        else:
            sdict[s]=None
    return mon(sdict)

def add_initial(model,pattern,value):
    import pysb.core
    complex_pattern = pysb.as_complex_pattern(pattern)

    for other_cp, other_value in model.initial_conditions:
        if complex_pattern.is_equivalent_to(other_cp):
            model.initial_conditions.remove((other_cp, other_value))
    model.initial(complex_pattern,value)
            

# Generate model rules via belpy
_,model = rdf_to_pysb.rdf_to_pysb('../../data/RAS_combined.rdf')

# Add ligand to the model
model.add_component(Monomer('EGF'))
model.add_component(Rule('EGF_activates_EGFR',model.monomers['EGF']()+model.monomers['EGFR']({'Kinase':'inactive'}) >> model.monomers['EGF']() + model.monomers['EGFR']({'Kinase':'active'}),model.parameters['kf_one_step_activate']))


# Add initial conditions
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

add_initial(model,model.monomers['EGF'](),model.parameters['EGF_0'])
add_initial(model,get_base_state(model.monomers['EGFR']),model.parameters['EGFR_0'])
add_initial(model,get_base_state(model.monomers['SOS1']),model.parameters['SOS_0'])
add_initial(model,get_base_state(model.monomers['HRAS']),model.parameters['RAS_0'])
add_initial(model,get_base_state(model.monomers['NRAS']),model.parameters['RAS_0'])
add_initial(model,get_base_state(model.monomers['KRAS']),model.parameters['RAS_0'])
add_initial(model,get_base_state(model.monomers['ARAF']),model.parameters['RAF_0'])
add_initial(model,get_base_state(model.monomers['BRAF']),model.parameters['RAF_0'])
add_initial(model,get_base_state(model.monomers['RAF1']),model.parameters['RAF_0'])
add_initial(model,get_base_state(model.monomers['MAP2K1']),model.parameters['MAP2K1_0'])
add_initial(model,get_base_state(model.monomers['MAPK1']),model.parameters['MAPK1_0'])
add_initial(model,model.monomers['RASA2'](Catalytic='active'),model.parameters['init_default'])
add_initial(model,model.monomers['RASA3'](Catalytic='active'),model.parameters['init_default'])

# Add observables
model.add_component(Observable("ERKact",model.monomers['MAPK1'](Kinase='active')))
model.add_component(Observable("MEKact",model.monomers['MAP2K1'](Kinase='active')))

# Import solver and generate model equations
#from pysb.integrate import Solver
#t = np.linspace(0,25,11)
#solver = Solver(model,t)
