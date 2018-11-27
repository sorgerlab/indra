import sys
import os
from os.path import join, dirname, abspath
from indra import preassembler
from indra.sources import eidos
from indra.sources.hume.make_hume_tsv import make_file as mht
from indra.sources.sofia.make_sofia_tsv import make_file as mst
from indra.java_vm import autoclass

eidos_package = 'org.clulab.wm.eidos'

if __name__ == '__main__':
    sofia_ont_path = sys.argv[1]
    hume_path = 'hume_ontology_examples.tsv'
    mht(hume_path)
    sofia_path = 'sofia_ontology_examples.tsv'
    mst(sofia_ont_path, sofia_path)

    om = autoclass(eidos_package + '.apps.OntologyMapper')
    eidos = autoclass(eidos_package + '.EidosSystem')
    es = eidos(autoclass('java.lang.Object')())

    example_weight = 0.8
    parent_weight = 0.1
    topn = 10
    table_str = om.mapOntologies(es, hume_path, sofia_path, example_weight,
                                 parent_weight, topn)
    with open(join(dirname(abspath(__file__)), os.pardir, 'resources',
                   'wm_ontomap.tsv'), 'w') as fh:
        fh.write(table_str)
