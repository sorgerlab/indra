import sys
import re
import rdflib
import os
import subprocess
import pkg_resources
from processor import BelProcessor
from . import ndex



def process_ndex_neighborhood(gene_names):
    bel_script = ndex.query_to_belscript(gene_names)
    with open('tmp.bel', 'wt') as fh:
        fh.write(bel_script)
    bel2rdf_path = pkg_resources.resource_filename('indra.bel', 'bel2rdf.rb')
    bel2rdf_cmd = "ruby %s --bel tmp.bel > tmp.rdf" % bel2rdf_path
    with open('tmp.rdf', 'wt') as fh:
        subprocess.call(bel2rdf_cmd.split(' '), stdout=fh, stderr=subprocess.STDOUT)
    with open('tmp.rdf', 'rt') as fh:
        rdf = fh.read()
    res = re.findall(r'_:([^ ]+)', rdf)
    for r in res:
        rdf = rdf.replace(r, r.replace('-', ''))
    with open('tmp2.rdf', 'w') as fh:
        fh.write(rdf)
    bp = process_belrdf('tmp2.rdf')
    bp.print_statements()
    return bp


def process_belrdf(rdf_filename):
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
    return bp

if __name__ == '__main__':
    # Make sure the user passed in an RDF filename
    if len(sys.argv) < 2:
        print "Usage: python bel_api.py file.rdf"
        sys.exit()
    # We take the RDF filename as the argument
    rdf_filename = sys.argv[1]
    bp = process_belrdf(rdf_filename)
