from jpype import *
import os
import sys

def getHGNC(proteinName):
    return 'http://purl.org/pc2/6/RelationshipXref_hgnc_symbol_%s_identity' % proteinName

def getSignature(e):
    if isinstance(e,JClass("org.biopax.paxtools.model.level3.BiochemicalReaction")):
        s = ' + '.join([getSignature(x) for x in e.left])
        if e.controlledOf:
            s += ' ==(' + '|'.join([getSignature(x) for x in e.controlledOf]) + ')==> '
        else:
            s += ' ====> '
        s += ' + '.join([getSignature(x) for x in e.right])
        return s
    elif isinstance(e,JClass("org.biopax.paxtools.model.level3.Catalysis")):
        s = ','.join([getSignature(x) for x in e.controller])
        return s
    elif isinstance(e,JClass("org.biopax.paxtools.model.level3.Protein")):
        s = e.displayName
        s += '(' + ','.join([f.toString() for f in e.feature]) + ')'
        return s
    elif isinstance(e,JClass("org.biopax.paxtools.model.level3.Complex")):
        s = '/'.join([getSignature(c) for c in e.component])
        return s
    else:
        if e.displayName:
            return e.displayName
        else:
            return ''


# Set default the location of the paxtools jar
# If CLASSPATH is set, this will not be used
# This jar was downloaded from
# http://sourceforge.net/projects/biopax/files/paxtools/
default_paxtools_jar = 'paxtools-4.3.0.jar'
# This example data file has been extracted and renamed from
# http://www.pathwaycommons.org/pc2/downloads/Pathway%20Commons.6.NCI%20Pathway%20Interaction%20Database:%20Pathway.BIOPAX.owl.gz
data_file = '../data/pathwaycommons_nci.owl'

JVM_memuse = '10g'

# Start the JVM
java_classpath = os.environ.get('CLASSPATH')
if java_classpath is None:
    print "CLASSPATH environment variable not set, using default classpath %s" % default_paxtools_jar
    java_classpath = default_paxtools_jar
startJVM(getDefaultJVMPath(), "-ea", "-Xmx%s" % JVM_memuse, "-Djava.class.path=%s" % java_classpath)

#get the paxtools root package as a shortcut
paxPkg = JPackage("org.biopax.paxtools")
io = paxPkg.io.SimpleIOHandler(paxPkg.model.BioPAXLevel.L3)

#import a BioPAX model from data_file
try:
   fileIS = java.io.FileInputStream(data_file)
   model = io.convertFromOWL(fileIS)
   fileIS.close()
except JavaException:
    print 'Could not read data file %s' % data_file
    sys.exit(0)

# These are the proteins whose neighborhood we are interested in
query_proteins = ['BRAF',]
#query_proteins = ['HRAS','KRAS','NRAS','ARAF','BRAF','RAF1','MAP2K1','MAPK1']

print 'Querying model for %s' % ','.join(query_proteins)

# Construct a set of the BioPax model elements corresponding to the proteins
query_set = JPackage('java.util').HashSet()
for p in query_proteins:
    pe = model.getByID(getHGNC(p))
    if pe is not None:
        query_set.add(pe)
    else:
        print 'Could not find protein %s in model' % p

filters = JArray(paxPkg.query.wrapperL3.Filter, 1)(0)
# Execute neighborhood query
result_set = paxPkg.query.QueryExecuter.runNeighborhood(query_set,model,1,paxPkg.query.algorithm.Direction.BOTHSTREAM,filters)

print 'Neighborhood query returned %d elements' % result_set.size()
print '-----------------'

for r in result_set:
    if isinstance(r,JClass("org.biopax.paxtools.model.level3.BiochemicalReaction")):
        print getSignature(r)

shutdownJVM() 
