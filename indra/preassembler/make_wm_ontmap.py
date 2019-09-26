import os
import sys
import requests
from os.path import join, dirname, abspath
from indra.sources.sofia.make_sofia_tsv import make_file as mst
from indra.java_vm import autoclass

eidos_package = 'org.clulab.wm.eidos'
wm_metadata_url = ('https://raw.githubusercontent.com/'
                   'WorldModelers/Ontologies/master/wm_metadata.yml')

if __name__ == '__main__':
    sofia_ont_path = sys.argv[1]
    sofia_path = 'sofia_ontology_examples.tsv'
    mst(sofia_ont_path, sofia_path)

    om = autoclass(eidos_package + '.apps.OntologyMapper')
    eidos = autoclass(eidos_package + '.EidosSystem')
    es = eidos()

    example_weight = 0.8
    parent_weight = 0.1
    topn = 5
    wm_yml = requests.get(wm_metadata_url).text
    ont_arg = autoclass('scala.Some')(wm_yml)
    table_str = om.mapOntologies(es, sofia_path, ont_arg, 'WM',
                                 example_weight, parent_weight, topn)
    with open(join(dirname(abspath(__file__)), os.pardir, 'resources',
                   'wm_ontomap.tsv'), 'w') as fh:
        fh.write(table_str)
