import rdflib
from os import pardir
from os.path import abspath, dirname, join
from rdflib import Namespace, Literal
from indra.databases.obo_client import OboClient


rdf_file = join(dirname(abspath(__file__)), pardir, 'resources',
                'cellular_component_hierarchy.rdf')


def save_hierarchy(g, path):
    with open(path, 'wb') as out_file:
        g_bytes = g.serialize(format='nt')
        # Replace extra new lines in string and get rid of empty line at end
        g_bytes = g_bytes.replace(b'\n\n', b'\n').strip()
        # Split into rows and sort
        rows = g_bytes.split(b'\n')
        rows.sort()
        g_bytes = b'\n'.join(rows)
        out_file.write(g_bytes)


def make_component_hierarchy(obo_client):
    g = rdflib.Graph()
    indra_ns = 'http://sorger.med.harvard.edu/indra/'
    #ln = Namespace('https://identifiers.org/')
    rn = Namespace(indra_ns + 'relations/')
    ln = Namespace(indra_ns + 'locations/')
    part_of = rn.term('partof')
    has_name = rn.term('hasName')
    for go_id, entry in obo_client.entries.items():
        if not entry['namespace'] == 'cellular_component':
            continue
        g.add((ln.term(go_id), has_name, Literal(entry['name'])))
        for parent_id in entry['is_a']:
            g.add((ln.term(go_id), part_of, ln.term(parent_id)))
    return g


if __name__ == '__main__':
    cl = OboClient(prefix='go')
    g = make_component_hierarchy(cl)
    save_hierarchy(g, rdf_file)
