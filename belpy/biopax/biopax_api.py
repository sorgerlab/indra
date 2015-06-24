import sys
import re
import subprocess

from processor import BiopaxProcessor

from jnius import autoclass

def owl_to_model(fname):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3) 
    
    try:
        fileIS = autoclass('java.io.FileInputStream')(fname)
    except JavaException:
        print 'Could not open data file %s' % fname
        sys.exit(0)
    try:
        biopax_model = io.convertFromOWL(fileIS)
    except JavaException:
        print 'Could not convert data file %s to BioPax model' % data_file
        sys.exit(0)
    
    fileIS.close()

    return biopax_model

def process_pc_neighborhood(gene_names):
    bp.print_statements()
    cpath_client = autoclass('cpath.client.CPathClient').newInstance('http://www.pathwaycommons.org/pc2/')
    query = cpath_client.createGraphQuery()
    query.kind(autoclass('cpath.service.GraphType').PATHSBETWEEN)
    query.sources(query_proteins)
    query.organismFilter(['homo sapiens'])
    query.mergeEquivalentInteractions(True)
    query.limit(autoclass('java.lang.Integer')(2))
    # Execute query
    model = query.result()
    bproc = BiopaxProcessor(model)
    bproc.get_complexes()
    bproc.print_statements()
    return bproc

def process_owl(owl_filename):
    model = owl_to_model(owl_filename)
    bproc = BiopaxProcessor(model)
    bproc.get_complexes()
    bproc.print_statements()
    return bproc

if __name__ == '__main__':
    # Make sure the user passed in an OWL filename
    if len(sys.argv) < 2:
        print "Usage: python biopax_api.py file.owl"
        sys.exit()
    # We take the OWL filename as the argument
    owl_filename = sys.argv[1]
    bproc = process_owl(owl_filename)
