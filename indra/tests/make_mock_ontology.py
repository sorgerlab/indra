"""This script can be used to create a mock bio ontology
which can be put in the appropriate cache location in place of the
real bio ontology for testing purposes"""
import os
import pickle
from indra.ontology.bio.ontology import BioOntology, CACHE_DIR


always_include = {
    'FPLX:ERK', 'HGNC:6871', 'HGNC:6877',
    'FPLX:AKT', 'FPLX:RAF', 'FPLX:MEK', 'FPLX:AMPK',
    'FPLX:SHC', 'FPLX:MAPK', 'FPLX:JNK',
    'FPLX:FOS_family', 'FPLX:JUN_family',
    'HGNC:9376', 'HGNC:9377', 'HGNC:9378', 'HGNC:9379',
    'HGNC:9385', 'HGNC:9386', 'HGNC:9387', 'FPLX:SRC', 'HGNC:391', 'HGNC:9955',
    'HGNC:6840', 'HGNC:6871', 'UP:Q13422',
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
    'UP:Q8BGM7', 'EFO:0000694', 'GO:GO:0005783', 'UP:Q13422',
    'MESH:D000938', 'FPLX:HIF_alpha', 'FPLX:HIF',
    'CHEBI:CHEBI:87307', 'CHEBI:CHEBI:36962',
    'UP:P15056', 'UP:Q32ZE1', 'UP:P15056', 'UP:P28482', 'UP:Q6P5R6',
    'UP:P62993', 'HGNC:4566', 'HGNC:18181', 'HGNC:10840',
    'HGNC:29869', 'HGNC:16743', 'GO:GO:0005737', 'GO:GO:0005575',
    'GO:GO:0005622', 'CHEBI:CHEBI:22950', 'CHEBI:CHEBI:37581',
    'CHEBI:CHEBI:25000', 'CHEBI:CHEBI:35701', 'CHEBI:CHEBI:36963',
    'GO:GO:0005886', 'GO:GO:0005737', 'GO:GO:0098826',
    'GO:GO:0016020', 'GO:GO:0005634',
    'UP:Q02750', 'UP:P01112', 'UP:P01019', 'UP:Q9MZT7', 'UP:Q13422',
    'HMDB:HMDB0000122', 'HGNC:7', 'HGNC:5', 'MIRBASE:MI0001730',
    'HGNC:31476', 'DRUGBANK:DB00001', 'MESH:D013812', 'CHEBI:CHEBI:26523',
    'UP:Q99490', 'MESH:D008099'
}

always_include_ns = {'FPLX', 'INDRA_ACTIVITIES', 'INDRA_MODS'}


def keep_node(node):
    ns = bio_ontology.get_ns(node)
    if ns in always_include_ns:
        return True
    if node in always_include:
        return True
    neigh = set(bio_ontology.successors(node)) | \
        set(bio_ontology.predecessors(node))
    if neigh & always_include:
        return True
    if {bio_ontology.get_ns(n) for n in neigh} & always_include_ns:
        return True
    return False


if __name__ == '__main__':
    bio_ontology = BioOntology()
    bio_ontology.initialize()
    keep_nodes = set()
    for node in bio_ontology.nodes:
        ns = bio_ontology.get_ns(node)
        if keep_node(node):
            keep_nodes.add(node)
    for node in list(bio_ontology.nodes):
        if node not in keep_nodes:
            bio_ontology.remove_node(node)
    bio_ontology._build_name_lookup()
    bio_ontology._label_components()
    bio_ontology._build_transitive_closure()
    fname = os.path.join(CACHE_DIR, 'mock_ontology.pkl')
    with open(fname, 'wb') as fh:
        pickle.dump(bio_ontology, fh, protocol=pickle.HIGHEST_PROTOCOL)
