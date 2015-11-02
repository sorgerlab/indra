#!/usr/bin/env python

import requests
import pandas as pd
import json
from jsonschema import validate
from . import ndexSchema as schema

        
# Utility to strip UUID from property graph before writing, ensure that we create new network
def removeUUIDFromNetwork(ndexPropertyGraphNetwork):   
    counter = 0
    for pv in ndexPropertyGraphNetwork.properties:
        if pv.predicateString == "UUID":
            del ndexPropertyGraphNetwork.properties[counter]
            return None
        else:
            counter = counter + 1
   
# Each ndex property becomes a dict property
def addNdexPropertiesToDict(properties, target):
     for pv in properties:
         predicate = pv['predicateString']  
         #dataType = pv['dataType']
         value = pv['value']
         # next, convert values based on datatype...
         target[predicate] = value
  
class PDGraph:
    def __init__(self, ndexPropertyGraphNetwork):
        # assemble a dictionary of nodes and then create the dataframe 
        node_rows = []
        node_indexes = []
        for index, node in ndexPropertyGraphNetwork['nodes'].iteritems():
            node_dict = {}
            # PropertyGraphNodes just have properties
            addNdexPropertiesToDict(node['properties'], node_dict)
            node_rows.append(node_dict)
            node_indexes.append(index)
            
        self.nodes = pd.DataFrame(node_rows, node_indexes)
        
        # assemble an array of dictionaries of edges where the edge id is the key
        # and the dictionary properties are the edge properties, the predicate, and the 
        # ids for subject and object.
        # then create the dataframe
        edge_rows = []
        edge_indexes = []
        for index, edge in ndexPropertyGraphNetwork['edges'].iteritems():
            edge_dict = {}
            # PropertyGraphEdges have special fields describing the 
            # edge itself, in addition to properties of the edge
            edge_dict['ndex:subjectId'] = edge['subjectId']
            edge_dict['ndex:predicate'] = edge['predicate']
            edge_dict['ndex:objectId'] = edge['objectId']
            addNdexPropertiesToDict(edge['properties'], edge_dict)
            edge_rows.append(edge_dict)
            edge_indexes.append(index)
            
        self.edges = pd.DataFrame(edge_rows, edge_indexes)
        
        # assemble a dictionary of the network properties
        property_dict = {}
        addNdexPropertiesToDict(ndexPropertyGraphNetwork['properties'], property_dict)
        self.properties = pd.Series(property_dict)
        
        ## TBD export a PDGraph to an NDEx PropertyGraphNetwork dict structure

class Ndex:
        
    def __init__(self, host = "http://www.ndexbio.org", username = None, password = None):
        if "localhost" in host:
            self.host = "http://localhost:8080/ndexbio-rest"
        else:
            self.host = host + "/rest"
        # create a session for this Ndex
        self.s = requests.session()
        if username and password:
            # add credentials to sesson, if available
            self.s.auth = (username, password)
    
# Base methods for making requests to this NDEx
    
    def put(self, route, putJson):
        url = self.host + route  
        print "PUT route: " + url
        print putJson
        headers = {'Content-Type' : 'application/json;charset=UTF-8',
                   'Accept' : 'application/json',
                   'Cache-Control' : 'no-cache',
                   }
        response = self.s.put(url, data = putJson, headers = headers)
        response.raise_for_status()
        return response.json()
        
    def post(self, route, postJson):
        url = self.host + route
        print "POST route: " + url
        print postJson
        
        headers = {'Content-Type': 'application/json',
                   'Accept': 'application/json',
                   'Cache-Control': 'no-cache',
                   }
        response = self.s.post(url, data=postJson, headers=headers)
        response.raise_for_status()
        return response.json()
        
    def delete(self, route):
        url = self.host + route
        response = self.s.delete(url)
        response.raise_for_status()
        return response.json()
    
    def get(self, route, getParams = None):
        url = self.host + route
        print "GET route: " + url
        print getParams
        response = self.s.get(url, params = getParams)
        response.raise_for_status()
        return response.json()
        
# Network methods
        

# Search for networks by keywords
#    network    POST    /network/search/{skipBlocks}/{blockSize}    SimpleNetworkQuery    NetworkSummary[]
    def findNetworks(self, searchString="", accountName=None, skipBlocks=0, blockSize=100): 
        route = "/network/search/%s/%s" % (skipBlocks, blockSize)
        postData = {"searchString" : searchString}
        if accountName:
            postData["accountName"] = accountName
        postJson = json.dumps(postData)
        return self.post(route, postJson)
        
    def findNetworksAsDataFrame(self, searchString="", accountName=None, skipBlocks=0, blockSize=100): 
        return pd.DataFrame(self.findNetworks(searchString, accountName, skipBlocks, blockSize))

    def getNetworkApi(self, asType="DataFrame"):
        route = "/network/api"
        decodedJson = self.get(route)
        if asType == "DataFrame":
            return pd.DataFrame(decodedJson)
        else:
            return decodedJson
 
#    network    POST    /network/{networkUUID}/edge/asNetwork/{skipBlocks}/{blockSize}        Network
    def getNetworkByEdges(self, networkId, skipBlocks=0, blockSize=100):
        route = "/network/%s/edge/asNetwork/%s/%s" % (networkId, skipBlocks, blockSize)
        return self.get(route)

#    network    GET    /network/{networkUUID}/asNetwork       Network
    def getCompleteNetwork(self, networkId):
        route = "/network/%s/asNetwork" % (networkId)
        return self.get(route)

#    network    GET    /network/{networkUUID}       NetworkSummary
    def getNetworkSummary(self, networkId):
        route = "/network/%s" % (networkId)
        return self.get(route)
        
#    network    POST    /network    Network    NetworkSummary
    def saveNewNetwork(self, Network):
        route = "/network/asNetwork"
        return self.post(route, Network)

#    network    POST    /network/asNetwork/group/{group UUID}    Network    NetworkSummary
    def saveNewNetworkForGroup(self, Network, groupId):
        route = "/network/asNetwork/group/%s" % (groupId)
        self.removeUUIDFromNetwork(Network)
        return self.post(route, Network)
       
##  Neighborhood PathQuery
#    network    POST    /network/{networkUUID}/asPropertyGraph/query    SimplePathQuery    PropertyGraphNetwork    
    def getNeighborhood(self, networkId, searchString, searchDepth=1):
        route = "/network/%s/asNetwork/query" % (networkId) 
        postData = {'searchString': searchString,
                   'searchDepth': searchDepth}
        postJson = json.dumps(postData)
        return self.post(route, postJson)
        
# PropertyGraphNetwork methods
        
#    network    POST    /network/{networkUUID}/edge/asPropertyGraph/{skipBlocks}/{blockSize}        PropertyGraphNetwork
    def getPropertyGraphNetworkByEdges(self, networkId, skipBlocks=0, blockSize=100):
        route = "/network/%s/edge/asPropertyGraph/%s/%s" % (networkId, skipBlocks, blockSize)
        return self.get(route)

#    network    GET    /network/{networkUUID}/asPropertyGraph        PropertyGraphNetwork
    def getCompletePropertyGraphNetwork(self, networkId):
        route = "/network/%s/asPropertyGraph" % (networkId)
        return self.get(route)

#    network    POST    /network/asPropertyGraph    PropertyGraphNetwork    NetworkSummary
    def saveNewPropertyGraphNetwork(self, propertyGraphNetwork):
        route = "/network/asPropertyGraph"
        self.removeUUIDFromNetwork(propertyGraphNetwork)
        return self.post(route, propertyGraphNetwork)

#    network    POST    /network/asPropertyGraph/group/{group UUID}    PropertyGraphNetwork    NetworkSummary
    def saveNewPropertyGraphNetworkForGroup(self, propertyGraphNetwork, groupId):
        route = "/network/asPropertyGraph/group/%s" % (groupId)
        self.removeUUIDFromNetwork(propertyGraphNetwork)
        return self.post(route, propertyGraphNetwork)
       
##  Neighborhood PathQuery
#    network    POST    /network/{networkUUID}/asPropertyGraph/query    SimplePathQuery    PropertyGraphNetwork    
    def getNeighborhoodAsPropertyGraph(self, networkId, searchString, searchDepth=1):
        route = "/network/%s/asPropertyGraph/query" % (networkId) 
        postData = {'searchString': searchString,
                   'searchDepth': searchDepth}
        postJson = json.dumps(postData)
        return self.post(route, postJson)

        
# User methods
        
# Group methods

# Request methods

# Task methods

# Validation

    def validateNetwork(self, network):
        return validate(network, schema.network)
