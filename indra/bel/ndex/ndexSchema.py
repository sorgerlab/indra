# -*- coding: utf-8 -*-
"""
Created on Sat Oct 18 17:04:44 2014

@author: Dexter Pratt
"""

network = {
    "$schema": "http://www.ndexbio.org/api/schema/1.0.0/Network",
    "description": "Schema for Network objects used in NDEx REST API",
    "type": "object",
    "required": [
        "type", 
        "edges", 
        "nodes", 
        "baseTerms", 
        "functionTerms", 
        "reifiedEdgeTerms", 
        "citations", 
        "supports", 
        "properties", 
        "presentationProperties",
        "namespaces", 
        "name"
        ],
    "properties": {
        "type": {
            "description": "NDEx object type. Must always = 'Network' for Network objects",
            "type": "string"
        },
        "creationTime": {
            "description": "Timestamp indicating Networks creation time. Ignored and can be null when Network is POSTed for creation.",
            "type": ["number", "null"],
        },
        "modificationTime": {
            "description": "Timestamp indicating Networks last modification time. Ignored and can be null when Network is POSTed or PUT",
            "type": ["number", "null"],
        },
        "isComplete": {
            "description": "Set to false while the network is being incrementally created or modified, true otherwise. Ignored and can be null when Network is POSTed or PUT",
            "type": ["boolean", "null"],
        },
        "isLocked": {
            "description": "Content modification permitted only if false. Ignored and can be null when Network is POSTed or PUT",
            "type": ["boolean", "null"],
        },
        "externalId": {
            "description": "UUID for the Network. Ignored and can be null when Network is POSTed for creation",
            "type": ["string", "null"]
        },
        "uri": {
            "description": "Unique URI for the Network based on the server address and Network UUID. Ignored and can be null when Network is POSTed for creation",
            "type": "string"
        },
        "nodeCount": {
            "description": "Number of Node objects in the Network, set by server when Network is returned. Ignored and can be null when Network is POSTed or PUT",
            "type": ["integer", "null"]
        },
        "edgeCount": {
            "description": "Number of Edge objects in the Network, set by server when Network is returned. Ignored and can be null when Network is POSTed or PUT",
            "type": "integer"
        },
        "visibility": {
            "description": "Visibility status of the Network. Set by server when Network is returned. Ignored and can be null when Network is POSTed or PUT",
            "type": ["string", "null"],
            "enum": [
                "PUBLIC",
                "PRIVATE",
                "DISCOVERABLE"
            ]
        },
        "name": {
            "description": "Name or Title of the Network. Required, but not a unique identifier for the Network.",
            "type": "string"
        },
        "description": {
            "description": "Text description of the Network.",
            "type": ["string", "null"]
        },
        "version": {
            "description": "Version string for the Network. No required format but strings conforming to semantic versioning standard are recommended and may be parsed.",
            "type": ["string", "null"]
        },
        "namespaces": {
            "description": "Object keys are integer element ids, values are Namespace objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/Namespace"
            },   
        },
        "baseTerms": {
            "description": "Object keys are integer element ids, values are BaseTerm objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/BaseTerm"
            }   
        },
        "functionTerms": {
            "description": "Object keys are integer element ids, values are FunctionTerm objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/FunctionTerm"
            }
        },
        "reifiedEdgeTerms": {
            "description": "Object keys are integer element ids, values are ReifiedEdgeTerm objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/ReifiedEdgeTerm"
            }
        },
        "nodes": {
            "description": "Object keys are integer element ids, values are Node objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/Node"
            }
        },
        "edges": {
            "description": "Object keys are integer element ids, values are Edge objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/Edge"
            }
        },
        "citations": {
            "description": "Object keys are integer element ids, values are Citation objects",
            "type": "object",
            "items": {
                "$ref": "#/definitions/Citation"
            }
        },
        "supports": {
            "description": "Object keys are integer element ids, values are Support objects",
            "type": "object",
            "items":{
                "$ref": "#/definitions/Support"
            }
        },
        "presentationProperties": {
            "description": "Items are SimplePropertyValuePair objects describing the appearance and layout of the Network",
            "type": "array",
            "items": {
                "$ref": "#/definitions/SimplePropertyValuePair"
            }
        },
        "properties": {
            "description": "Items are NdexPropertyValuePair objects describing the content of the Network",
            "type": "array",
            "items": {
                "$ref": "#/definitions/NdexPropertyValuePair"            
            }
        }
    },
    "definitions": {
        "Edge": {
            "type": "object",
            "required": [
                "type", 
                "id", 
                "citationIds", 
                "supportIds",
                "subjectId",
                "properties", 
                "presentationProperties"
                ],
            "properties": {
                "type": {
                    "type": "string"
                },
                "id": {
                    "type": "integer"
                },
                "subjectId": {
                    "type": "integer"
                },
                "predicateId": {
                    "type": ["integer", "null"]
                },
                "objectId": {
                    "type": ["integer", "null"]
                },
                "presentationProperties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/SimplePropertyValuePair"
                    }
                },
                "properties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/NdexPropertyValuePair"           
                    }
                },
                "citationIds": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }
                },
                "supportIds": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }
                }
            }   
        },
        "Node": {
            "type": "object",
            "required": [
                "type", 
                "id", 
                "citationIds", 
                "supportIds",
                "aliases",
                "relatedTerms",
                "properties", 
                "presentationProperties"
                ],
            "properties": {
                "type": {
                    "type": "string"
                },
                "id": {
                    "type": "integer"
                },
                "name": {
                    "type": ["string", "null"]
                },
                "aliases": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }       
                },
                "citationIds": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }
                },
                "supportIds": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }
                },
                "presentationProperties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/SimplePropertyValuePair"
                    }
                },
                "properties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/NdexPropertyValuePair"             
                    }
                },
                "relatedTerms": {
                    "type": "array",
                    "items": {
                        "type": "integer"
                    }
                },
                "represents": {
                    "type": ["integer", "null"]
                },
                "representsTermType": {
                    "type": ["string", "null"]
                }
            }          
        },
        "Namespace": {
            "type" : "object",
            "required": [
                "type", 
                "id", 
                "properties", 
                "presentationProperties"
                ],
            "properties" : {
                "type" : {
                    "type" : "string"
                },
                "id" : {
                    "type" : "integer"
                },
                "prefix" : {
                    "type" : ["string", "null"]
                },
                "uri" : {
                    "type" : ["string", "null"]
                },
                "presentationProperties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/SimplePropertyValuePair"
                    }
                },
                "properties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/NdexPropertyValuePair"              
                    }
                }
            }
        },
        "BaseTerm": {
            "type" : "object",
            "required": [
                "type", 
                "id",
                "name",
                "termType"
                ],
            "properties" : {
                "id" : {
                    "type" : "integer"
                },
                "type" : {
                    "type" : "string"
                },
                "name" : {
                    "type" : "string"
                },
                "namespaceId" : {
                    "type" : ["integer", "null"]
                },
                "termType" : {
                    "type" : "string"
                }
            }
        },
        "FunctionTerm": {
            "type" : "object",
            "required": [
                "type", 
                "id", 
                "termType", 
                "functionTermId",
                "parameterIds"
                ],
            "properties" : {
                "id" : {
                    "type" : "integer"
                },
                "type" : {
                    "type" : "string"
                },
                "parameterIds" : {
                    "type" : "array",
                    "items" : {
                        "type" : "integer"
                    }
                },
                "functionTermId" : {
                    "type" : "integer"
                },
                "termType" : {
                    "type" : "string"
                }
            }
        },
        "ReifiedEdgeTerm": {
            "type" : "object",
            "required": [
                "type", 
                "id", 
                "termType", 
                "presentationProperties"
                ],
            "properties" : {
                "id" : {
                    "type" : "integer"
                },
                "type" : {
                    "type" : "string"
                },
                "termType" : {
                    "type" : "string"
                },
                "edgeId" : {
                    "type" : "integer"
                }
            }
        },
        "Citation": {
            "type" : "object",
            "required": [
                "type", 
                "id",
                "contributors",
                "properties", 
                "presentationProperties",
                "identifier",
                "idType"
                ],
            "properties" : {
                "id" : {
                    "type" : "integer"
                },
                "type" : {
                    "type" : "string"
                },
                "contributors" : {
                    "type" : "array",
                    "items" : {
                        "type" : "string"
                    }
                },
                "title" : {
                    "type" : ["string", "null"]
                },
                "identifier" : {
                    "type" : "string"
                },
                "idType" : {
                    "type" : "string"
                },
                "presentationProperties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/SimplePropertyValuePair"
                    }
                },
                "properties": {
                    "type": "array",
                    "items": {
                        "$ref": "#/definitions/NdexPropertyValuePair"              
                    }
                }
            }                
        },
        "Support": {
            "type" : "object",
            "required": [
                "type", 
                "id", 
                "text"
                ],
            "properties" : {
                "id" : {
                    "type" : "integer"
                },
                "type" : {
                    "type" : "string"
                },
                "text" : {
                    "type" : "string"
                },
                "citationId" : {
                    "type" : ["integer", "null"]
                }
            }
        },
        "SimplePropertyValuePair":{
            "type" : "object",
            "required": [
                "type", 
                "name", 
                "value"
                ],
            "properties" : {
                "type" : {
                    "type" : "string"
                },
                "name" : {
                    "type" : "string"
                },
                "value" : {
                    "type" : "string"
                }
            }
        },
        "NdexPropertyValuePair":{
            "type" : "object",
            "required": [
                "type"
                ],
            "properties" : {
                "type" : {
                    "type" : "string"
                },
                "value" : {
                    "type" : ["string", "null"]
                },
                "valueId" : {
                    "type" : ["integer", "null"]
                },
                "dataType" : {
                    "type" : ["string", "null"]
                },
                "predicateString" : {
                    "type" : ["string", "null"]
                },
                "predicateId" : {
                    "type" : ["integer", "null"]
                }
            }
        }
    }
}