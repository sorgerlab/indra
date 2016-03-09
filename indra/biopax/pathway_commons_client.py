import urllib, urllib2
from indra.java_vm import autoclass, JavaException

pc2_url = 'http://www.pathwaycommons.org/pc2/'

def send_request(kind, source, target=None):
    kind_str = kind.lower()
    if kind not in ['neighborhood', 'pathsbetween', 'pathsfromto']:
        print 'Invalid query type %s' % kind_str
        return None
    organism = '9606'
    if isinstance(source, basestring):
        source_str = source
    else:
        source_str = ','.join(source)
    params = {'kind': kind_str,
              'organism': organism,
              'source': ','.join(source),
              'format': 'BIOPAX'}
    if target is not None:
        if isinstance(target, basestring):
            target_str = target
        else:
            target_str = ','.join(target)
        params['target'] = target_str
    
    print 'Sending Pathway Commons query...'
    res = urllib2.urlopen(pc2_url + 'graph', data=urllib.urlencode(params))
    owl_str = res.read()
    model = owl_str_to_model(owl_str)
    if model is not None:
        print 'Pathway Commons query returned a model...'
    return model

def owl_str_to_model(owl_str):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)
    bais = autoclass('java.io.ByteArrayInputStream')
    scs = autoclass('java.nio.charset.StandardCharsets')
    jstr = autoclass('java.lang.String')
    istream = bais(jstr(owl_str).getBytes(scs.UTF_8));
    biopax_model = io.convertFromOWL(istream)
    return biopax_model

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
