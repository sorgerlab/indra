import os
import sys
import ipdb

# The CLASSPATH environmental variable needs to be set
# to point to the paxtools and the cpath jar files.
# These can be obtained from the sites below
# http://sourceforge.net/projects/biopax/files/paxtools/
# https://code.google.com/p/pathway-commons/wiki/PC2Client
import jnius_config
jnius_config.add_options('-Xmx4g')
from jnius import autoclass, JavaException

use_data_file = False

def getHGNC(proteinName):
    return 'http://purl.org/pc2/7/RelationshipXref_hgnc_symbol_%s_identity' % proteinName

def getSignature(e):
    if isinstance(e,autoclass("org.biopax.paxtools.impl.level3.BiochemicalReactionImpl")):
        s = '\nFrom: '
        s += ' + '.join([getSignature(x) for x in e.getLeft().toArray()])
        if e.getControlledOf():
            catalyzers = [x for x in e.getControlledOf().toArray()
                if isinstance(x,autoclass("org.biopax.paxtools.impl.level3.CatalysisImpl"))]
            controllers = [x for x in e.getControlledOf().toArray() 
                if isinstance(x,autoclass("org.biopax.paxtools.impl.level3.ControlImpl"))]
            if catalyzers:
                s += '\nCatalyzed by: '
                s += ' | '.join([getSignature(x) for x in catalyzers])
            if controllers:
                s += '\nControlled by: '
                s += ' | '.join([getSignature(x) for x in controllers])
        s += '\nTo: ' + ' + '.join([getSignature(x) for x in e.getRight().toArray()])
        xrefs = e.getXref().toArray()
        if xrefs:
            s += '\nReferences: '
            s += ' | '.join([x.getRDFId() for x in xrefs])
        data_source = e.getDataSource().toArray()[0]
        if data_source:
            s += '\nData source: '
            s += data_source.getDisplayName()
            data_source_comment = data_source.getComment()
            if data_source_comment:
                s += '('+data_source_comment.toArray()[0]+')'
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.EvidenceImpl")):
        s = ' | '.join(e.getComment().toArray())
        s += ', ' + ' | '.join([getSignature(x) for x in e.getEvidenceCode().toArray()])
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.EvidenceCodeVocabularyImpl")):
        return e.toString()
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.CatalysisImpl")):
        s = ','.join([getSignature(x) for x in e.getController().toArray()])
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ProteinImpl")):
        s = e.getDisplayName()
        s += '(' + ';'.join([getSignature(f) for f in e.getFeature().toArray()]) + ')'
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ComplexImpl")):
        #s = '/'.join([getSignature(c) for c in e.getComponent().toArray()])
        phys_ent = e.getMemberPhysicalEntity().toArray()
        compnt = e.getComponent().toArray()
        s = ''
        if phys_ent:
            s += '%'.join([getSignature(c) for c in phys_ent])
        elif compnt:
            s += '%'.join([getSignature(c) for c in compnt])
        else:
            s += e.getDisplayName()
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ModificationFeatureImpl")):
        term_acc = autoclass("org.biopax.paxtools.controller.PathAccessor")("ModificationFeature/modificationType/term")
        site_acc = autoclass("org.biopax.paxtools.controller.PathAccessor")("ModificationFeature/featureLocation:SequenceSite/sequencePosition")
        feat_acc = autoclass("org.biopax.paxtools.controller.PathAccessor")("PhysicalEntity/feature:ModificationFeature")
        
        terms = ','.join(term_acc.getValueFromBean(e).toArray())
        sites = ','.join(['%d'%x for x in site_acc.getValueFromBean(e).toArray()])
        feats = feat_acc.getValueFromBean(e).toArray()
        s = terms
        if sites != '':
            s += '@' + sites
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.EntityFeatureImpl")):
        mods = [c if not c.startswith('REPLACED') else c for c in e.getComment().toArray()]
        s = ';'.join(mods)
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ControlImpl")):
        s = ','.join([getSignature(x) for x in e.getController().toArray()])
        return s
    else:
        if hasattr(e,'getDisplayName') and e.getDisplayName():
            return e.getDisplayName()
        else:
            return e.toString()

def saveResult(result,file_name):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)
    try:
        fileOS = autoclass('java.io.FileOutputStream')(file_name)
        io.convertToOWL(result, fileOS)
        fileOS.close()
    except JavaException:
        print 'Could not write to file %s' % file_name

# These are the proteins whose neighborhood we are interested in
query_proteins = ['EGFR','SHC1','GRB2','SOS1','NRAS','HRAS','KRAS','RAF1','ARAF','BRAF','MAP2K1','MAP2K2']

if use_data_file:
    # This example data file has been extracted and renamed from
    # http://www.pathwaycommons.org/archives/PC2/v6-201502/Pathway%20Commons.6.NCI%20Pathway%20Interaction%20Database:%20Pathway.BIOPAX.owl.gz
    #data_file = '../data/pathwaycommons_nci.owl'
    data_file = '/home/beni/data/pathwaycommons/pathwaycommons_detailed.owl'
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    print 'Starting offline query of %s' % data_file
    #import a BioPAX model from data_file
    try:
        fileIS = autoclass('java.io.FileInputStream')(data_file)
        model = io.convertFromOWL(fileIS)
        fileIS.close()
    except JavaException:
        print 'Could not read data file %s' % data_file
        sys.exit(0)
    # Construct a set of the BioPax model elements corresponding to the proteins
    query_set = autoclass('java.util.HashSet')()
    for p in query_proteins:
        pe = model.getByID(getHGNC(p))
        if pe is not None:
            query_set.add(pe)
        else:
            print 'Could not find protein %s in model' % p

    filters = autoclass('org.biopax.paxtools.query.wrapperL3.Filter')
    qe = autoclass('org.biopax.paxtools.query.QueryExecuter')
    # Execute query
    #result_set = qe.runNeighborhood(query_set,model,1,autoclass('org.biopax.paxtools.query.algorithm.Direction').BOTHSTREAM)
    result_set = qe.runPathsBetween(query_set,model,1)
else:
    print 'Starting online query of the Pathway Commons database'
    cpath_client = autoclass('cpath.client.CPathClient').newInstance('http://www.pathwaycommons.org/pc2/')
    query = cpath_client.createGraphQuery()
    query.kind(autoclass('cpath.service.GraphType').PATHSBETWEEN)
    query.sources(query_proteins)
    query.organismFilter(['homo sapiens'])
    query.mergeEquivalentInteractions(True)
    query.limit(autoclass('java.lang.Integer')(2))
    # Execute query
    result_model = query.result()
    result_set = result_model.getObjects()
    saveResult(result_model,'test.owl')

print 'Queried model for %s' % ','.join(query_proteins)

reactions = [r for r in result_set.toArray() if 
    isinstance(r,autoclass("org.biopax.paxtools.impl.level3.BiochemicalReactionImpl"))]

print 'Query returned %d elements (%d reactions)' % (result_set.size(),len(reactions))
print '-----------------'


# Pretty print all reactions
for r in reactions:
    print getSignature(r)
        
