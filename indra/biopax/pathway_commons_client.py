from indra.java_vm import autoclass, JavaException

def run_pc_query(query_type, source_genes, target_genes=None, neighbor_limit=1):
    cpath_client = autoclass('cpath.client.CPathClient').\
        newInstance('http://www.pathwaycommons.org/pc2/')
    query = cpath_client.createGraphQuery()
    query.kind(query_type)
    query.sources(source_genes)
    query.targets(target_genes)
    query.organismFilter(['homo sapiens'])
    query.mergeEquivalentInteractions(True)
    query.limit(autoclass('java.lang.Integer')(neighbor_limit))
    # Execute query
    print 'Sending Pathway Commons query...'
    model = query.result()
    if model is not None:
        print 'Pathway Commons query returned model...'
    else:
        print 'Pathway Commons query returned blank model...'
    return model

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
