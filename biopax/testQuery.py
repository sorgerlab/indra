import os
import sys
import ipdb

# The CLASSPATH environmental variable needs to be set
# to point to the paxtools and the cpath jar files.
# These can be obtained from the sites below
# http://sourceforge.net/projects/biopax/files/paxtools/
# https://code.google.com/p/pathway-commons/wiki/PC2Client
import jnius_config
jnius_config.add_options('-Xmx3000m')
from jnius import autoclass, JavaException

use_data_file = True

def getHGNC(proteinName):
    return 'http://purl.org/pc2/6/RelationshipXref_hgnc_symbol_%s_identity' % proteinName

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
        #evidence =  e.getEvidence().toArray()
        #if evidence:
        #    s += '\nEvidence: '
        #    s += ' | '.join([getSignature(x) for x in evidence])
        #comments = e.getComment().toArray()
        #if comments:
        #    s += '\nComments: '
        #    s += ' | '.join(comments)
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
        s += '(' + ','.join([getSignature(f) for f in e.getFeature().toArray()]) + ')'
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ComplexImpl")):
        s = '/'.join([getSignature(c) for c in e.getComponent().toArray()])
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.EntityFeatureImpl")):
        s = ','.join(e.getComment().toArray())
        return s
    elif isinstance(e,autoclass("org.biopax.paxtools.impl.level3.ControlImpl")):
        s = ','.join([getSignature(x) for x in e.getController().toArray()])
        return s
    else:
        if hasattr(e,'getDisplayName') and e.getDisplayName():
            return e.getDisplayName()
        else:
            return e.toString()


# These are the proteins whose neighborhood we are interested in
query_proteins = ['BRAF']
#query_proteins = ['HRAS','KRAS','NRAS','ARAF','BRAF','RAF1','MAP2K1','MAPK1']


if use_data_file:
    # This example data file has been extracted and renamed from
    # http://www.pathwaycommons.org/archives/PC2/v6-201502/Pathway%20Commons.6.NCI%20Pathway%20Interaction%20Database:%20Pathway.BIOPAX.owl.gz
    data_file = '../data/pathwaycommons_nci.owl'
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)


    print 'Starting offline query of the NCI dataset'
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
    result_set = qe.runNeighborhood(query_set,model,1,autoclass('org.biopax.paxtools.query.algorithm.Direction').BOTHSTREAM)
else:
    print 'Starting online query of the Pathway Commons database'
    cp = autoclass('cpath.client.CPathClient').newInstance('http://www.pathwaycommons.org/pc2/')
    query = cp.createGraphQuery()
    query.sources(query_proteins)
    # Execute query
    result_model = query.result()
    result_set = result_model.getObjects()


print 'Querying model for %s' % ','.join(query_proteins)

print 'Neighborhood query returned %d elements' % result_set.size()
print '-----------------'

for r in result_set.toArray():
    if isinstance(r,autoclass("org.biopax.paxtools.impl.level3.BiochemicalReactionImpl")):
        print getSignature(r)
        
