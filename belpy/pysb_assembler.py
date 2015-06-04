from pysb import *
from pysb.core import SelfExporter
from bel import bel_api

SelfExporter.do_export = False

def add_default_initial_conditions(model):
    try:
        default_ic = model.parameters['default_ic']
    except KeyError:
        default_ic = Parameter('default_ic', 100.)
        model.add_component(default_ic)
    # Iterate over all monomers
    for m in model.monomers:
        # Build up monomer pattern dict
        sites_dict = {}
        for site in m.sites:
            if site in m.site_states:
                sites_dict[site] = m.site_states[site][0]
            else:
                sites_dict[site] = None
        mp = m(**sites_dict)
        model.initial(mp, default_ic)

class pysb_assembler(object):
    def __init__(self):
        self.statements = []
    
    def add_statements(self,stmts):
        self.statements.extend(stmts)

    def make_model(self,initial_conditions=True):
        model = Model()
        for stmt in self.statements:
            stmt.monomers(model)
        for stmt in self.statements:
            stmt.assemble(model)
        if initial_conditions:
            add_default_initial_conditions(model)
        return model

if __name__ == '__main__':
  pa = pysb_assembler()
  bp = bel_api.process_belrdf('data/RAS_neighborhood.rdf')
  pa.add_statements(bp.belpy_stmts)
  model = pa.make_model()
