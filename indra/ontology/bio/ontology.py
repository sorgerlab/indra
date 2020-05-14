import os
from ..ontology_graph import IndraOntology, label
from indra.util import read_unicode_csv
from indra.databases import hgnc_client, uniprot_client, chebi_client, \
    mesh_client, obo_client
from indra.sources.trips.processor import ncit_map
from indra.statements import modtype_conditions


HERE = os.path.dirname(os.path.abspath(__file__))
resources = os.path.join(HERE, os.pardir, os.pardir, 'resources')


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
        self.add_mesh_xrefs()
        # Add hierarchies
        self.add_famplex_hierarchy()
        self.add_obo_hierarchies()
        self.add_mesh_hierarchy()
        self.add_activity_hierarchy()
        self.add_modification_hierarchy()

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
        rel_mappings = {
            'xref': 'xref',
            'isa': 'isa',
            'partof': 'partof',
            'is_a': 'isa',
            'part_of': 'partof',
            # These are for ChEBI: identical to the old behavior but it might
            # make sense to add other relations here too
            'is_conjugate_acid_of': 'isa',
            'has_functional_parent': 'isa',
            'has_parent_hydride': 'isa',
            'has_role': 'isa'
        }
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                for rel, targets in entry.get('relations', {}).items():
                    # Skip unknown relation types
                    mapped_rel = rel_mappings.get(rel)
                    if not mapped_rel:
                        continue
                    for target in targets:
                        edges.append((label(ns.upper(), db_id),
                                      label(ns.upper(), target),
                                      {'type': mapped_rel}))
        self.add_edges_from(edges)

    def add_chemical_xrefs(self):
        mappings = [
            (chebi_client.chebi_chembl, 'CHEBI', 'CHEMBL'),
            (chebi_client.chebi_pubchem, 'CHEBI', 'PUBCHEM'),
            (chebi_client.hmdb_chebi, 'HMDB', 'CHEBI'),
            (chebi_client.cas_chebi, 'CAS', 'CHEBI'),
        ]
        edges = []
        data = {'type': 'xref', 'source': 'chebi'}

        def label_fix(ns, id):
            if ns == 'CHEBI' and not id.startswith('CHEBI'):
                id = 'CHEBI:%s' % id
            return label(ns, id)

        for map_dict, from_ns, to_ns in mappings:
            for from_id, to_id in map_dict.items():
                source = label_fix(from_ns, from_id)
                target = label_fix(to_ns, to_id)
                edges.append((source, target, data))
                edges.append((target, source, data))
        self.add_edges_from(edges)

    def add_mesh_nodes(self):
        nodes = [(label('MESH', mesh_id),
                  {'name': name})
                 for mesh_id, name in
                 mesh_client.mesh_id_to_name.items()]
        self.add_nodes_from(nodes)

    def add_mesh_xrefs(self):
        edges = []
        data = {'type': 'xref', 'source': 'gilda'}
        for mesh_id, (db_ns, db_id) in mesh_client.mesh_to_db.items():
            edges.append((label('MESH', mesh_id),
                          label(db_ns, db_id),
                          data))
        for (db_ns, db_id), mesh_id in mesh_client.db_to_mesh.items():
            edges.append((label(db_ns, db_id),
                          label('MESH', mesh_id),
                          data))
        self.add_edges_from(edges)

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

    def add_activity_hierarchy(self):
        rels = [
            ('transcription', 'activity'),
            ('catalytic', 'activity'),
            ('gtpbound', 'activity'),
            ('kinase', 'catalytic'),
            ('phosphatase', 'catalytic'),
            ('gef', 'catalytic'),
            ('gap', 'catalytic')
        ]
        self.add_edges_from([
            (label('INDRA_ACTIVITIES', source),
             label('INDRA_ACTIVITIES', target),
             {'type': 'isa'})
            for source, target in rels
        ]
        )

    def add_modification_hierarchy(self):
        self.add_edges_from([
            (label('INDRA_MODS', source),
             label('INDRA_MODS', 'modification'),
             {'type': 'isa'})
            for source in modtype_conditions
            if source != 'modification'
        ]
        )


bio_ontology = BioOntology()
