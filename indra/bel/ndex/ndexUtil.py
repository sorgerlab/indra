# -*- coding: utf-8 -*-
"""
Created on Sun Oct  5 11:10:59 2014

@author: Dexter Pratt
"""
import sys
import networkx as nx
import re
import StringIO

# Convert NDEx property graph json to a trivial networkx network
def ndexPropertyGraphNetworkToNetworkX(ndexPropertyGraphNetwork):
        g = nx.MultiDiGraph()
        for node in ndexPropertyGraphNetwork['nodes'].values():
            g.add_node(node['id'])
        for edge in ndexPropertyGraphNetwork['edges'].values():
            g.add_edge(edge['subjectId'], edge['objectId'])
        return g

# This is specific to summarizing a BEL network. 
# Need to generalize
def stripPrefixes(input):
    st = input.lower()
    if st.startswith('bel:'):
        return input[4:len(input)]
    elif st.startswith('hgnc:'):
         return input[5:len(input)]
    else:
         return input

# This is BEL specific, since BEL is the only current user of funciton terms
def getFunctionAbbreviation(input):
    st = input.lower()
    fl = stripPrefixes(st)
    if fl == "abundance":
        return "a"
    elif fl == "biological_process":
        return "bp"
    elif fl ==  "catalytic_activity":
        return "cat"
    elif fl ==  "complex_abundance":
        return "complex"
    elif fl ==  "pathology":
        return "path"
    elif fl ==  "peptidase_activity":
        return "pep"
    elif fl ==  "protein_abundance":
        return "p"
    elif fl ==  "rna_abundance":
        return "r"
    elif fl ==  "protein_modification":
        return "pmod"
    elif fl ==  "transcriptional_activity":
        return "tscript"
    elif fl ==  "molecular_activity":
        return "act"
    elif fl ==  "degradation":
        return "deg"
    elif fl ==  "kinase_activity":
        return "kin"
    elif fl ==  "substitution":
        return "sub"
    else:
        return fl

def getFunctionFull(input):
    st = input.lower()
    fl = stripPrefixes(st)
    if fl == "abundance":
        return "abundance"
    elif fl == "biologicalprocess":
        return "biologicalProcess"
    elif fl ==  "catalyticactivity":
        return "catalyticActivity"
    elif fl ==  "complexabundance":
        return "complexAbundance"
    elif fl ==  "pathology":
        return "pathology"
    elif fl ==  "peptidaseactivity":
        return "peptidaseActivity"
    elif fl ==  "proteinabundance":
        return "proteinAbundance"
    elif fl ==  "rnaabundance":
        return "rnaAbundance"
    elif fl ==  "proteinmodification":
        return "proteinModification"
    elif fl ==  "transcriptionalactivity":
        return "transcriptionalActivity"
    elif fl ==  "molecularactivity":
        return "molecularActivity"
    elif fl ==  "degradation":
        return "degradation"
    elif fl ==  "phosphataseactivity":
        return "phosphataseActivity"
    elif fl ==  "kinaseactivity":
        return "kinaseActivity"
    elif fl ==  "cellsecretion":
        return "cellSecretion"
    elif fl ==  "substitution":
        return "substitution"
    elif fl == "gtpboundactivity":
        return "gtpBoundActivity"
    elif fl == "cellsurfaceexpression":
        return "cellSurfaceExpression"
    elif fl == "micrornaabundance":
        fl = "microRnaAbundance"
    else:
        return fl

def getPredicateFull(p):
    if p == 'INCREASES':
        return 'increases'
    elif p == 'DECREASES':
        return 'decreases'
    elif p == 'DIRECTLY_INCREASES':
        return 'directlyIncreases'
    elif p == 'DIRECTLY_DECREASES':
        return 'directlyDecreases'
    else:
        return p

class NetworkWrapper:
    def __init__(self, ndexNetwork, removeNamespace=None):
        self.network = ndexNetwork
        self.supportToEdgeMap = {}
        self.citationToSupportMap = {}
        self.nodeLabelMap = {}
        self.termLabelMap = {}

        for nodeId, node in ndexNetwork['nodes'].iteritems():
            nodeLabel = self.getNodeLabel(node)
            removeNode = False
            for rn in removeNamespace:
                if nodeLabel.find(rn)!=-1:
                    removeNode = True

            if removeNode==False:
                self.nodeLabelMap[int(nodeId)] = self.getNodeLabel(node)

        for edge in ndexNetwork['edges'].values():
            for supportId in edge['supportIds']:
                supports = ndexNetwork['supports']
                support = supports[str(supportId)]
                if supportId in self.supportToEdgeMap:
                    edgeList = self.supportToEdgeMap[supportId]
                else:
                    edgeList = []
                edgeList.append(edge)
                self.supportToEdgeMap[supportId] = edgeList

        for supportId in self.supportToEdgeMap.keys():
            support = ndexNetwork['supports'][str(supportId)]
            citationId = support['citationId']
            if citationId in self.citationToSupportMap:
                supportIdList = self.citationToSupportMap[citationId]
            else:
                supportIdList = []
            supportIdList.append(supportId)
            self.citationToSupportMap[citationId] = supportIdList

    def getEdgeLabel(self, edge):
        subjectLabel = "missing"
        objectLabel = "missing"
        predicateLabel = "missing"
        subjectId = edge['subjectId']
        objectId = edge['objectId']
        if subjectId in self.nodeLabelMap:
            subjectLabel = self.nodeLabelMap[subjectId]
        if objectId in self.nodeLabelMap:
            objectLabel = self.nodeLabelMap[objectId]
        predicateId = edge['predicateId']
        predicateLabel = stripPrefixes(self.getTermLabel(predicateId))
        predicateLabel = getPredicateFull(predicateLabel)
        label = "%s %s %s" % (subjectLabel, predicateLabel, objectLabel)
        return label

    def getNodeLabel(self, node):
        if 'name' in node and node['name']:
            return node['name']

        elif 'represents' in node:
            return self.getTermLabel(node['represents'])

        else:
            return "node %s" % (node['id'])

    def getTermById(self, termId):
        termIdStr = str(termId)
        if termIdStr in self.network['baseTerms']:
            return self.network['baseTerms'][termIdStr]
        elif termIdStr in self.network['functionTerms']:
            return self.network['functionTerms'][termIdStr]
        elif termIdStr in self.network['reifiedEdgeTerms']:
            return self.network['reifiedEdgeTerms'][termIdStr]
        else:
            return None

    def getTermLabel(self, termId):
        if termId in self.termLabelMap:
            return self.termLabelMap[termId]
        else:
            label = "error"
            term = self.getTermById(termId)
            type = term['type'].lower()
            if type == "baseterm":
                name = term['name']
                if 'namespaceId' in term and term['namespaceId']:
                    namespaceId = term['namespaceId']
                    #namespace = self.network['namespaces'][namespaceId]
                    try:
                        namespace = self.network['namespaces'][str(namespaceId)]
                    except KeyError:
                        namespace = None

                    if namespace:
                        if namespace['prefix']:
                            if namespace['prefix']!='BEL':
                                if re.search('[^a-zA-Z0-9]',name) != None:
                                    name = '"'+name+'"'
                                label = "%s:%s" % (namespace['prefix'], name)
                            else:
                                label = name
                        elif namespace['uri']:
                            label = "%s%s" % (namespace['uri'], name)
                        else:
                            label = name
                    else:
                        label = name
                else:
                    label = name

            elif type == "functionterm":
                functionTermId = term['functionTermId']
                functionLabel = self.getTermLabel(functionTermId)
                #functionLabel = getFunctionAbbreviation(functionLabel)
                functionLabel = getFunctionFull(functionLabel)
                parameterLabels = []
                for parameterId in term['parameterIds']:
                    parameterLabel = self.getTermLabel(parameterId)
                    parameterLabels.append(parameterLabel)
                label = "%s(%s)" % (functionLabel, ",".join(parameterLabels))

            elif type == "reifiededgeterm":
                edgeId = term['edgeId']
                edges = self.network['edges']
                if edgeId in edges:
                    reifiedEdge = edges[edgeId]
                    label = "(%s)" % (self.getEdgeLabel(reifiedEdge))
                else:
                    label = "(reifiedEdge: %s)" % (edgeId)

            else:
                label = "term: %s" % (termId)

            self.termLabelMap[termId] = label
            return label

    def write_utf(self,output,string):
        output.write(string.encode('utf8','replace'))
        #output.write(string)
        
    def writeBELScript(self, fileName = None):
        output = StringIO.StringIO()
       
        self.write_utf(output,'#Properties section\n')
        self.write_utf(output,'SET DOCUMENT Name = "NDEx query result in BEL script"\n')
        self.write_utf(output,'SET DOCUMENT Description = "Query with ndex-python-client, one step neighborhood"\n')

        
        # Print definitions in header
        self.write_utf(output,'\n# Definitions Section\n')

        # Print namespaces
        for _,ns in self.network['namespaces'].iteritems():
            if ns['uri'].endswith('.belns'):
                self.write_utf(output,'DEFINE NAMESPACE %s AS URL "%s"\n' % (ns['prefix'],ns['uri']))
        
        # Print annotations
        for _,ann in self.network['namespaces'].iteritems():
            if ann['uri'].endswith('.belanno'):
                self.write_utf(output,'DEFINE ANNOTATION %s AS URL "%s"\n' % (ann['prefix'],ann['uri']))

        # Print BEL statements
        self.write_utf(output,'\n#Statements section\n')

        print 
        print 'Unhandled statements'
        print '===================='
    
        # Iterate by citation
        for citationId, supportIdList in self.citationToSupportMap.iteritems():
            # Start a group for each citation
            self.write_utf(output,'\nSET STATEMENT_GROUP = "Group %d"\n' % citationId)
            try:
                citation = self.network['citations'][str(citationId)]
                citation_title = citation['title']
                citation_terms = citation['identifier'].split(':')
                if citation_terms[0]=='pmid':
                    citation_type = 'PubMed'
                    citation_id = citation_terms[1]
                else:
                    citation_type = 'N/A'
                    citation_id = citation['identifier']
                self.write_utf(output,('SET Citation = {"%s","%s","%s"}\n' % (citation_type, citation_title, citation_id)))    
            except KeyError:
                self.write_utf(output,'SET Citation = {"","",""}\n')
           
            # Iterate by evidence within each citation
            for supportId in supportIdList:
                support = self.network['supports'][str(supportId)]
                supportText = support['text'].replace('"','').replace('\n',' ')
                self.write_utf(output,('\nSET Evidence = "%s"\n' % supportText))
                edgeList = self.supportToEdgeMap[supportId]
                # Print BEL statements 
                for edge in edgeList:
                    outstr = self.getEdgeLabel(edge)
                    if outstr.find('missing') != -1:
                        continue

                    # Generate valid translocation statements - not used
                    #outstr = re.sub(r'GOCCACC:GO:(\d+),GOCCACC:GO:(\d+)',r'fromLoc(GOCCACC:\1),toLoc(GOCCACC:\2)',outstr)

                    # Reified edges not handled
                    if outstr.find('reifiedEdge') == -1:
                        # Translocation not handled
                        if outstr.find('translocation') == -1:
                            # 'None' modifiers not handled
                            if outstr.find('None') == -1:
                                print outstr
                                self.write_utf(output,"%s\n" % outstr)
                            else:
                                print outstr
                        else:
                            print outstr
                    else:
                        print outstr
            self.write_utf(output,'\nUNSET STATEMENT_GROUP\n')

        retstr = output.getvalue()
        if fileName:
            outfile = open(fileName, 'wt')
            outfile.write(retstr)
            outfile.close()

        output.close()
        return retstr

    def writeSummary(self, fileName = None):
        if fileName:
            output = open(fileName, 'w')
        else:
            output = sys.stdout
            
        for citationId, supportIdList in self.citationToSupportMap.iteritems():
            citations = self.network['citations']
            citation = citations[str(citationId)]
            citationId = citation['identifier']
            # Write Citation
            output.write("\n=========================================================================\n")
            output.write("        Citation: %s\n" % (citationId))
            output.write("=========================================================================\n\n")

            for supportId in supportIdList:
                support = self.network['supports'][str(supportId)]
                # Write Support
                output.write("_______________________________\n")
                output.write(("Evidence: %s\n\n" % support['text']).encode('utf8','replace'))

                edgeList = self.supportToEdgeMap[supportId]
                for edge in edgeList:
                    # Write Edge
                    output.write("       %s\n" % self.getEdgeLabel(edge))
                    for pv in edge['properties']:
                        output.write("                %s: %s\n" % (pv['predicateString'], pv['value']))

        if fileName:
            output.close()


