from indra.sources import eidos
from indra.sources.hume.make_hume_tsv import make_file
from indra.java_vm import autoclass

eidos_package = 'org.clulab.wm.eidos'

if __name__ == '__main__':
    bbn_path = 'hume_examaples.tsv'
    make_file(bbn_path)
    sofia_path = 'sofia_examples.tsv'

    om = autoclass(eidos_package + '.apps.OntologyMapper')
    eidos = autoclass(eidos_package + '.EidosSystem')
    es = eidos(autoclass('java.lang.Object')())

    example_weight = 0.8
    parent_weight = 0.1
    topn = 10
    table_str = om.mapOntologies(es, bbn_path, sofia_path, example_weight,
                                 parent_weight, topn)
