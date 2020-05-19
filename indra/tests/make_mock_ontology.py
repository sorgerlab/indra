"""This script can be used to create a mock bio ontology
which oen can put in the appropriate cache location in place of the
real bio ontology for testing purposes"""
import os
import pickle
from collections import defaultdict
from indra.ontology.bio.ontology import bio_ontology, CACHE_DIR


always_include = {
    'FPLX:ERK', 'HGNC:6871', 'HGNC:6877',
    'FPLX:AKT', 'FPLX:RAF', 'FPLX:MEK', 'FPLX:AMPK',
    'HGNC:9376', 'HGNC:9377', 'HGNC:9378', 'HGNC:9379',
    'HGNC:9385', 'HGNC:9386', 'HGNC:9387', 'FPLX:SRC', 'HGNC:391', 'HGNC:9955',
    'CHEBI:CHEBI:76971', 'CHEBI:CHEBI:37045', 'CHEBI:CHEBI:15996',
    'CHEBI:CHEBI:75771', 'CHEBI:CHEBI:37121', 'CHEBI:CHEBI:57600',
    'UP:P04585', 'HP:HP:0031801', 'GO:GO:0006915',
    'MESH:D008545', 'MESH:D058750', 'GO:GO:0001837', 'MESH:D058750',
    'CHEBI:CHEBI:46661', 'MESH:D000067777', 'HGNC:3313', 'UP:Q12926',
    'HP:HP:0000002', 'DOID:DOID:0014667', 'EFO:1002050', 'EFO:0000001',
    'EFO:0009502', 'HP:HP:0031801', 'MESH:D064706', 'DOID:DOID:0060495',
    'MESH:D000071017', 'HP:HP:0031801', 'UPPRO:PRO_0000032458', 'HGNC:3467',
    'HGNC:13006', 'HGNC:6407', 'UP:Q15208', 'UP:Q92597', 'UP:Q6IE75',
    'CHEBI:CHEBI:63637', 'UP:P04608', 'UP:O43687', 'HGNC:377', 'UP:Q9UGI9',
    'UP:Q8BGM7', 'EFO:0000694', 'GO:GO:0005783',
}


if __name__ == '__main__':
    keep_nodes = defaultdict(list)
    for node in bio_ontology.nodes:
        if node in always_include or \
                set(bio_ontology.successors(node)) & always_include or \
                set(bio_ontology.predecessors(node)) & always_include:
            keep_nodes[ns].append(node)
        ns = bio_ontology.get_ns(node)
        if ns not in keep_nodes or len(keep_nodes.get(ns)) < 100:
            keep_nodes[ns].append(node)
    for node in list(bio_ontology.nodes):
        if node not in keep_nodes[bio_ontology.get_ns(node)]:
            bio_ontology.remove_node(node)
    fname = os.path.join(CACHE_DIR, 'mock_ontology.pkl')
    with open(fname, 'wb') as fh:
        pickle.dump(bio_ontology, fh, protocol=pickle.HIGHEST_PROTOCOL)
