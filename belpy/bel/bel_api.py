import sys
import re
import rdflib
from BelProcessor import BelProcessor
import subprocess

def process_ndex_neighborhood(gene_names):
    ndex = __import__('ndex-python-client')
    bel_script = ndex.query_to_belscript(gene_names)
    fh = open('tmp.bel', 'wt')
    fh.write(bel_script)
    fh.close()
    bel_to_rdf_cmd = "bel2rdf --bel tmp.bel > tmp.rdf"
    subprocess.call(bel_to_rdf_cmd.split(' '))
    fh = open('tmp.rdf', 'rt')
    rdf = fh.read()
    fh.close()
    res = re.findall(r'_:([^ ]+)', rdf)
    for r in res:
        rdf = rdf.replace(r,r.replace('-',''))
    fh = open('tmp2.rdf','wt')
    fh.write(rdf)
    fh.close()
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
        print "Usage: python rdf_to_pysb.py file.rdf"
        sys.exit()
    # We take the RDF filename as the argument
    rdf_filename = sys.argv[1]
    bp = process_belrdf(rdf_filename)
