import sys

from processor import BiopaxProcessor

from indra.java_vm import autoclass, JavaException

def owl_to_model(fname):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    try:
        file_is = autoclass('java.io.FileInputStream')(fname)
    except JavaException:
        print 'Could not open data file %s' % fname
        return
    try:
        biopax_model = io.convertFromOWL(file_is)
    except JavaException:
        print 'Could not convert data file %s to BioPax model' % data_file
        return

    file_is.close()

    return biopax_model

def model_to_owl(model, fname):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    try:
        fileOS = autoclass('java.io.FileOutputStream')(fname)
    except JavaException:
        print 'Could not open data file %s' % fname
        return
    l3_factory = autoclass('org.biopax.paxtools.model.BioPAXLevel').L3.getDefaultFactory()
    model_out = l3_factory.createModel()
    for r in model.getObjects().toArray():
        model_out.add(r)
    io.convertToOWL(model_out, fileOS)

    fileOS.close()



def process_pc_neighborhood(gene_names, neighbor_limit=1):
    query_type = autoclass('cpath.service.GraphType').NEIGHBORHOOD
    model = _run_pc_query(query_type, gene_names, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_pc_pathsbetween(gene_names, neighbor_limit=1):
    query_type = autoclass('cpath.service.GraphType').PATHSBETWEEN
    model = _run_pc_query(query_type, gene_names, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_owl(owl_filename):
    model = owl_to_model(owl_filename)
    return process_model(model)


def process_model(model):
    bproc = BiopaxProcessor(model)
    # bproc.get_complexes()
    # bproc.get_phosphorylation()
    # bproc.print_statements()
    return bproc


def _run_pc_query(query_type, gene_names, neighbor_limit=1):
    cpath_client = autoclass('cpath.client.CPathClient').\
        newInstance('http://www.pathwaycommons.org/pc2/')
    query = cpath_client.createGraphQuery()
    query.kind(query_type)
    query.sources(gene_names)
    query.organismFilter(['homo sapiens'])
    query.mergeEquivalentInteractions(True)
    query.limit(autoclass('java.lang.Integer')(neighbor_limit))
    # Execute query
    model = query.result()
    return model

if __name__ == '__main__':
    # Make sure the user passed in an OWL filename
    if len(sys.argv) < 2:
        print "Usage: python biopax_api.py file.owl"
        sys.exit()
    # We take the OWL filename as the argument
    owl_filename = sys.argv[1]
    bproc = process_owl(owl_filename)
