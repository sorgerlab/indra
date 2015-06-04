import sys

import rdflib

from pysb import *
from pysb.core import SelfExporter

from belpy.BelProcessor import BelProcessor

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

def bel_rdf_to_pysb(rdf_filename, initial_conditions=True):
    # Parse the RDF
    g = rdflib.Graph()
    g.parse(rdf_filename, format='nt')
    # Build BelPy statements from RDF
    bp = BelProcessor(g)
    bp.get_complexes()
    bp.get_activating_subs()
    bp.get_modifications()
    bp.get_dephosphorylations()
    bp.get_activating_mods()
    bp.get_composite_activating_mods()
    bp.get_activity_activity()

    # Print some output about the process
    bp.print_statement_coverage()
    print "\n--- Converted BelPy Statements -------------"
    bp.print_statements()

    # Make the PySB model and return
    model = bp.make_model()
    if initial_conditions:
        add_default_initial_conditions(model)
    return (bp, model)

if __name__ == '__main__':
    # Make sure the user passed in an RDF filename
    if len(sys.argv) < 2:
        print "Usage: python rdf_to_pysb.py file.rdf"
        sys.exit()
    # We take the RDF filename as the argument
    rdf_filename = sys.argv[1]
    (bp, model) = bel_rdf_to_pysb(rdf_filename, initial_conditions=True)
