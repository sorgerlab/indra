import os
import networkx
from indra.util import read_unicode_csv
from indra.databases import hgnc_client, uniprot_client, chebi_client, \
    mesh_client, obo_client

HERE = os.path.dirname(os.path.abspath(__file__))


class IndraOntology(networkx.MultiDiGraph):
    def __init__(self):
        super().__init__()
        # Add all nodes with annotations
        self.add_hgnc_nodes()
        self.add_uniprot_nodes()
        self.add_famplex_nodes()
        self.add_obo_nodes()
        # Add xrefs
        self.add_hgnc_uniprot_xrefs()
        self.add_famplex_xrefs()
        # Add hierarchies
        self.add_famplex_hierarchy()
        self.add_obo_hierarchies()

    def add_hgnc_nodes(self):
        nodes = [(label('HGNC', hid),
                  {'name': hname, 'ns': 'HGNC', 'id': hid})
                 for (hid, hname) in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        nodes = [(label('UP', uid), {'name': uname, 'ns': 'UP', 'id': uid})
                 for (uid, uname)
                 in uniprot_client.um.uniprot_gene_name.items()]
        self.add_nodes_from(nodes)

    def add_hgnc_uniprot_xrefs(self):
        edges = []
        for hid, uid in hgnc_client.uniprot_ids.items():
            uids = uid.split(', ')
            for uid in uids:
                edges.append((label('HGNC', hid), label('UP', uid),
                              {'type': 'xref', 'source': 'hgnc'}))
        self.add_edges_from(edges)

        edges = [(label('UP', uid), label('HGNC', hid),
                  {'type': 'xref', 'source': 'hgnc'})
                 for uid, hid in uniprot_client.um.uniprot_hgnc.items()]
        self.add_edges_from(edges)

    def add_famplex_nodes(self):
        nodes = []
        for row in read_unicode_csv(os.path.join(HERE, 'famplex',
                                                 'entities.csv'),
                                    delimiter=','):
            entity = row[0]
            nodes.append((label('FPLX', entity),
                          {'name': entity, 'id': entity}))
        self.add_nodes_from(nodes)

    def add_famplex_hierarchy(self):
        edges = []
        for row in read_unicode_csv(os.path.join(HERE, 'famplex',
                                                 'relations.csv'),
                                    delimiter=','):
            ns1, id1, rel, ns2, id2 = row
            if ns1 == 'HGNC':
                id1 = hgnc_client.get_hgnc_id(id1)
            edges.append((label(ns1, id1), label(ns2, id2), {'type': rel}))
        self.add_edges_from(edges)

    def add_famplex_xrefs(self):
        edges = []
        include_refs = {'PF', 'IP', 'GO', 'NCIT'}
        for row in read_unicode_csv(os.path.join(HERE, 'famplex',
                                                 'equivalences.csv'),
                                    delimiter=','):
            ref_ns, ref_id, fplx_id = row
            if ref_ns not in include_refs:
                continue
            edges.append((label(ref_ns, ref_id),
                          label('FPLX', fplx_id),
                          {'type': 'xref', 'source': 'fplx'}))
            edges.append((label('FPLX', fplx_id),
                          label(ref_ns, ref_id),
                          {'type': 'xref', 'source': 'fplx'}))
        self.add_edges_from(edges)

    def add_obo_nodes(self):
        namespaces = ['go', 'efo', 'hp', 'efo']
        nodes = []
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                nodes.append((label(ns.upper(), db_id),
                              {'name': entry['name'], 'id': db_id}))
        self.add_nodes_from(nodes)

    def add_obo_hierarchies(self):
        namespaces = ['go', 'efo', 'hp', 'efo']
        edges = []
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                for rel, targets in entry.get('relations', {}).items():
                    for target in targets:
                        edges.append((label(ns.upper(), db_id),
                                      target, {'type': rel}))
        self.add_edges_from(edges)


def label(ns, id):
    return '%s:%s' % (ns, id)
