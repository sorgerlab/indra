import os
import logging
import networkx
from collections import deque
from indra.util import read_unicode_csv
from indra.databases import hgnc_client, uniprot_client, chebi_client, \
    mesh_client, obo_client
from indra.sources.trips.processor import ncit_map
from .shortest_path import bidirectional_shortest_path

HERE = os.path.dirname(os.path.abspath(__file__))
resources = os.path.join(HERE, os.pardir, 'resources')

logger = logging.getLogger(__name__)


class IndraOntology(networkx.MultiDiGraph):
    def _check_path(self, ns1, id1, ns2, id2, edge_types):
        try:
            path = bidirectional_shortest_path(self,
                                               label(ns1, id1),
                                               label(ns2, id2),
                                               edge_types=edge_types)
        except networkx.NetworkXError:
            return False
        return True

    @staticmethod
    def get_ns_id(node):
        return node.split(':', maxsplit=1)

    @staticmethod
    def get_ns(node):
        return IndraOntology.get_ns_id(node)[0]

    @staticmethod
    def get_id(node):
        return IndraOntology.get_ns_id(node)[1]

    def isrel(self, ns1, id1, ns2, id2, rels):
        return self._check_path(ns1, id1, ns2, id2, rels)

    def isa(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'isa'})

    def partof(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'partof'})

    def isa_or_partof(self, ns1, id1, ns2, id2):
        return self.isrel(ns1, id1, ns2, id2, rels={'isa', 'partof'})

    def maps_to(self, ns1, id1, ns2, id2):
        return self._check_path(ns1, id1, ns2, id2, {'xref'})

    def map_to(self, ns1, id1, ns2):
        source = label(ns1, id1)
        targets = [target for target in self.xrefs(source)
                   if self.get_ns(target) == ns2]
        if len(targets) == 1:
            return targets[0]
        return None

    def xrefs(self, source):
        for _, target, data in self.edges(source, data=True):
            if data['type'] == 'xref':
                yield target

    def get_mappings(self, ns, id):
        source = label(ns, id)
        visited = {source}
        queue = deque([(source, self.xrefs(source))])
        targets = []
        while queue:
            parent, children = queue[0]
            try:
                child = next(children)
                if child not in visited:
                    targets.append(child)
                    visited.add(child)
                    queue.append((child, self.xrefs(child)))
            except StopIteration:
                queue.popleft()
        return targets

    def get_name(self, ns, id):
        node = label(ns, id)
        try:
            return self.nodes[node]['name']
        except KeyError:
            return None


class BioOntology(IndraOntology):
    def __init__(self):
        super().__init__()
        # Add all nodes with annotations
        self.add_hgnc_nodes()
        self.add_uniprot_nodes()
        self.add_famplex_nodes()
        self.add_obo_nodes()
        self.add_mesh_nodes()
        self.add_ncit_nodes()
        # Add xrefs
        self.add_hgnc_uniprot_xrefs()
        self.add_famplex_xrefs()
        self.add_chemical_xrefs()
        self.add_ncit_xrefs()
        # Add hierarchies
        self.add_famplex_hierarchy()
        self.add_obo_hierarchies()
        self.add_mesh_hierarchy()


    def add_hgnc_nodes(self):
        nodes = [(label('HGNC', hid), {'name': hname})
                 for (hid, hname) in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        nodes = [(label('UP', uid), {'name': uname})
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
        for row in read_unicode_csv(os.path.join(resources, 'famplex',
                                                 'entities.csv'),
                                    delimiter=','):
            entity = row[0]
            nodes.append((label('FPLX', entity),
                          {'name': entity}))
        self.add_nodes_from(nodes)

    def add_famplex_hierarchy(self):
        edges = []
        for row in read_unicode_csv(os.path.join(resources, 'famplex',
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
        for row in read_unicode_csv(os.path.join(resources, 'famplex',
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
        namespaces = ['go', 'efo', 'hp', 'doid', 'chebi']
        nodes = []
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                nodes.append((label(ns.upper(), db_id),
                              {'name': entry['name']}))
        self.add_nodes_from(nodes)

    def add_obo_hierarchies(self):
        namespaces = ['go', 'efo', 'hp', 'doid', 'chebi']
        edges = []
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                for rel, targets in entry.get('relations', {}).items():
                    for target in targets:
                        edges.append((label(ns.upper(), db_id),
                                      label(ns.upper(), target),
                                      {'type': rel}))
        self.add_edges_from(edges)

    def add_chemical_xrefs(self):
        edges = []
        # Chebi/Chembl
        for chebi_id, chembl_id in chebi_client.chebi_chembl.items():
            edges.append((label('CHEBI', chebi_id),
                          label('CHEMBL', chembl_id),
                          {'type': 'xref', 'source': 'chebi'}))
            edges.append((label('CHEMBL', chembl_id),
                          label('CHEBI', chebi_id),
                          {'type': 'xref', 'source': 'chebi'}))

        # Chebi/PubChem
        for chebi_id, pubchem_id in chebi_client.chebi_pubchem.items():
            edges.append((label('CHEBI', chebi_id),
                          label('PUBCHEM', pubchem_id),
                          {'type': 'xref', 'source': 'chebi'}))
            edges.append((label('PUBCHEM', pubchem_id),
                          label('CHEBI', chebi_id),
                          {'type': 'xref', 'source': 'chebi'}))

        # Chebi/HMDB
        for hmdb_id, chebi_id in chebi_client.hmdb_chebi.items():
            edges.append((label('CHEBI', chebi_id),
                          label('HMDB', hmdb_id),
                          {'type': 'xref', 'source': 'chebi'}))
            edges.append((label('HMDB', hmdb_id),
                          label('CHEBI', chebi_id),
                          {'type': 'xref', 'source': 'chebi'}))

        # Chebi/CAS
        for cas_id, chebi_id in chebi_client.cas_chebi.items():
            edges.append((label('CHEBI', chebi_id),
                          label('CAS', cas_id),
                          {'type': 'xref', 'source': 'chebi'}))
            edges.append((label('CAS', cas_id),
                          label('CHEBI', chebi_id),
                          {'type': 'xref', 'source': 'chebi'}))

        self.add_edges_from(edges)

    def add_mesh_nodes(self):
        nodes = [(label('MESH', mesh_id),
                  {'name': name})
                 for mesh_id, name in
                 mesh_client.mesh_id_to_name.items()]
        self.add_nodes_from(nodes)

    def add_mesh_hierarchy(self):
        mesh_tree_numbers_to_id = {}
        for mesh_id, tns in mesh_client.mesh_id_to_tree_numbers.items():
            for tn in tns:
                mesh_tree_numbers_to_id[tn] = mesh_id
        edges = []
        for mesh_id, tns in mesh_client.mesh_id_to_tree_numbers.items():
            parents_added = set()
            for tn in tns:
                if '.' not in tn:
                    continue
                parent_tn, _ = tn.rsplit('.', maxsplit=1)
                parent_id = mesh_tree_numbers_to_id[parent_tn]
                if parent_id in parents_added:
                    continue
                edges.append((label('MESH', mesh_id),
                              label('MESH', parent_id),
                              {'type': 'isa'}))
        self.add_edges_from(edges)

    def add_ncit_nodes(self):
        nodes = [(label('NCIT', ncit_id)) for ncit_id in ncit_map]
        self.add_nodes_from(nodes)

    def add_ncit_xrefs(self):
        edges = []
        for ncit_id, (target_ns, target_id) in ncit_map.items():
            edges.append((label('NCIT', ncit_id),
                          label(target_ns, target_id),
                          {'type': 'xref', 'source': 'ncit'}))
        self.add_edges_from(edges)


def label(ns, id):
    return '%s:%s' % (ns, id)


bio_ontology = BioOntology()
