from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import requests
from indra.java_vm import autoclass, JavaException
# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('biopax')

pc2_url = 'http://www.pathwaycommons.org/pc2/'

def graph_query(kind, source, target=None, neighbor_limit=1):
    """Perform a graph query on PathwayCommons.

    For more information on these queries, see
    http://www.pathwaycommons.org/pc2/#graph

    Parameters
    ----------
    kind : str
        The kind of graph query to perform. Currently 3 options are
        implemented, 'neighborhood', 'pathsbetween' and 'pathsfromto'.
    source : list[str]
        A list of gene names which are the source set for the graph query.
    target : Optional[list[str]]
        A list of gene names which are the target set for the graph query.
        Only needed for 'pathsfromto' queries.
    neighbor_limit : Optional[int]
        This limits the length of the longest path considered in
        the graph query. Default: 1

    Returns
    -------
    model : org.biopax.paxtools.model.Model
        A BioPAX model (java object).
    """
    params = {}
    params['format'] = 'BIOPAX'
    params['organism'] = '9606'
    # Get the "kind" string
    kind_str = kind.lower()
    if kind not in ['neighborhood', 'pathsbetween', 'pathsfromto']:
        logger.warn('Invalid query type %s' % kind_str)
        return None
    params['kind'] = kind_str
    # Get the source string
    if isinstance(source, basestring):
        source_str = source
    else:
        source_str = ','.join(source)
    params['source'] = source_str
    try:
        neighbor_limit = int(neighbor_limit)
        params['limit'] = neighbor_limit
    except (TypeError, ValueError):
        logger.warn('Invalid neighborhood limit %s' % neighbor_limit)
        return None
    if target is not None:
        if isinstance(target, basestring):
            target_str = target
        else:
            target_str = ','.join(target)
        params['target'] = target_str

    logger.info('Sending Pathway Commons query...')
    res = requests.get(pc2_url + 'graph', params=params)
    if not res.status_code == 200:
        logger.error('Response is HTTP code %d.' % res.status_code)
        if res.status_code == 500:
            logger.error('Note: HTTP code 500 can mean empty '
                         'results for a valid query.')
        return None
    # We don't decode to Unicode here because owl_str_to_model expects
    # a byte stream
    model = owl_str_to_model(res.content)
    if model is not None:
        logger.info('Pathway Commons query returned a model...')
    return model

def owl_str_to_model(owl_str):
    """Return a BioPAX model object from an OWL string.

    Parameters
    ----------
    owl_str : str
        The model as an OWL string.

    Returns
    -------
    biopax_model : org.biopax.paxtools.model.Model
        A BioPAX model object (java object).
    """
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)
    bais = autoclass('java.io.ByteArrayInputStream')
    scs = autoclass('java.nio.charset.StandardCharsets')
    jstr = autoclass('java.lang.String')
    istream = bais(jstr(owl_str).getBytes(scs.UTF_8));
    biopax_model = io.convertFromOWL(istream)
    return biopax_model

def owl_to_model(fname):
    """Return a BioPAX model object from an OWL file.

    Parameters
    ----------
    fname : str
        The name of the OWL file containing the model.

    Returns
    -------
    biopax_model : org.biopax.paxtools.model.Model
        A BioPAX model object (java object).
    """
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    try:
        file_is = autoclass('java.io.FileInputStream')(fname)
    except JavaException:
        logger.error('Could not open data file %s' % fname)
        return
    try:
        biopax_model = io.convertFromOWL(file_is)
    except JavaException as e:
        logger.error('Could not convert data file %s to BioPax model' % fname)
        logger.error(e)
        return

    file_is.close()

    return biopax_model

def model_to_owl(model, fname):
    """Save a BioPAX model object as an OWL file.

    Parameters
    ----------
    model : org.biopax.paxtools.model.Model
        A BioPAX model object (java object).
    fname : str
        The name of the OWL file to save the model in.
    """
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    try:
        fileOS = autoclass('java.io.FileOutputStream')(fname)
    except JavaException:
        logger.error('Could not open data file %s' % fname)
        return
    l3_factory = autoclass('org.biopax.paxtools.model.BioPAXLevel').L3.getDefaultFactory()
    model_out = l3_factory.createModel()
    for r in model.getObjects().toArray():
        model_out.add(r)
    io.convertToOWL(model_out, fileOS)

    fileOS.close()
