"""This script can be used to create a mock bio ontology
which oen can put in the appropriate cache location in place of the
real bio ontology for testing purposes"""
import pickle
from collections import defaultdict
from indra.ontology.bio.ontology import bio_ontology


if __name__ == '__main__':
    keep_nodes = defaultdict(list)
    for node in bio_ontology.nodes:
        ns = bio_ontology.get_ns(node)
        if ns not in keep_nodes or len(keep_nodes.get(ns)) < 1000:
            keep_nodes[ns].append(node)
    for node in list(bio_ontology.nodes):
        if node not in keep_nodes[bio_ontology.get_ns(node)]:
            bio_ontology.remove_node(node)
    with open('mock_ontology.pkl', 'wb') as fh:
        pickle.dump(bio_ontology, fh, protocol=pickle.HIGHEST_PROTOCOL)
