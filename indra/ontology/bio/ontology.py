import os
import time
import pickle
import logging
from indra.config import get_config
from ..ontology_graph import IndraOntology
from indra.util import read_unicode_csv
from indra.databases import hgnc_client, uniprot_client, chebi_client, \
    mesh_client, obo_client
from indra.sources.trips.processor import ncit_map
from indra.statements import modtype_conditions


HERE = os.path.dirname(os.path.abspath(__file__))
resources = os.path.join(HERE, os.pardir, os.pardir, 'resources')
CACHE_DIR = get_config('INDRA_RESOURCES') or \
    os.path.join(os.path.expanduser('~'), '.indra')
CACHE_FILE = os.path.join(CACHE_DIR, 'bio_ontology.pkl')


logger = logging.getLogger(__name__)


class BioOntology(IndraOntology):
    def __init__(self):
        super().__init__()
        self._initialized = False

    def initialize(self):
        logger.info('Initializing bio ontology...')
        # Add all nodes with annotations
        ts = time.time()
        self.add_hgnc_nodes()
        self.add_uniprot_nodes()
        self.add_famplex_nodes()
        self.add_obo_nodes()
        self.add_mesh_nodes()
        self.add_ncit_nodes()
        self.add_uppro_nodes()
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
        self.add_uppro_hierarchy()
        self._initialized = True
        # Build name to ID lookup
        self._build_name_lookup()
        logger.info('Finished initializing bio ontology...')

    def add_hgnc_nodes(self):
        nodes = [(self.label('HGNC', hid), {'name': hname})
                 for (hid, hname) in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        nodes = [(self.label('UP', uid), {'name': uname})
                 for (uid, uname)
                 in uniprot_client.um.uniprot_gene_name.items()]
        self.add_nodes_from(nodes)

    def add_uppro_nodes(self):
        nodes = []
        for prot_id, features in uniprot_client.um.features.items():
            for feature in features:
                node = self.label('UPPRO', feature.id)
                data = {'name': feature.name}
                nodes.append((node, data))
        self.add_nodes_from(nodes)

    def add_hgnc_uniprot_xrefs(self):
        edges = []
        for hid, uid in hgnc_client.uniprot_ids.items():
            uids = uid.split(', ')
            for uid in uids:
                edges.append((self.label('HGNC', hid), self.label('UP', uid),
                              {'type': 'xref', 'source': 'hgnc'}))
        self.add_edges_from(edges)

        edges = [(self.label('UP', uid), self.label('HGNC', hid),
                  {'type': 'xref', 'source': 'hgnc'})
                 for uid, hid in uniprot_client.um.uniprot_hgnc.items()]
        self.add_edges_from(edges)

    def add_famplex_nodes(self):
        nodes = []
        for row in read_unicode_csv(os.path.join(resources, 'famplex',
                                                 'entities.csv'),
                                    delimiter=','):
            entity = row[0]
            nodes.append((self.label('FPLX', entity),
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
            edges.append((self.label(ns1, id1),
                          self.label(ns2, id2),
                          {'type': rel}))
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
            edges.append((self.label(ref_ns, ref_id),
                          self.label('FPLX', fplx_id),
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
                nodes.append((self.label(ns.upper(), db_id),
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
                        edges.append((self.label(ns.upper(), db_id),
                                      self.label(ns.upper(), target),
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
            return self.label(ns, id)

        for map_dict, from_ns, to_ns in mappings:
            for from_id, to_id in map_dict.items():
                source = label_fix(from_ns, from_id)
                target = label_fix(to_ns, to_id)
                edges.append((source, target, data))
                edges.append((target, source, data))
        self.add_edges_from(edges)

    def add_mesh_nodes(self):
        nodes = [(self.label('MESH', mesh_id),
                  {'name': name})
                 for mesh_id, name in
                 mesh_client.mesh_id_to_name.items()]
        self.add_nodes_from(nodes)

    def add_mesh_xrefs(self):
        edges = []
        data = {'type': 'xref', 'source': 'gilda'}
        for mesh_id, (db_ns, db_id) in mesh_client.mesh_to_db.items():
            edges.append((self.label('MESH', mesh_id),
                          self.label(db_ns, db_id),
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
                edges.append((self.label('MESH', mesh_id),
                              self.label('MESH', parent_id),
                              {'type': 'isa'}))
        self.add_edges_from(edges)

    def add_ncit_nodes(self):
        nodes = [(self.label('NCIT', ncit_id)) for ncit_id in ncit_map]
        self.add_nodes_from(nodes)

    def add_ncit_xrefs(self):
        edges = []
        for ncit_id, (target_ns, target_id) in ncit_map.items():
            edges.append((label('NCIT', ncit_id),
                          label(target_ns, target_id),
                          {'type': 'xref', 'source': 'ncit'}))
        self.add_edges_from(edges)

    def add_uppro_hierarchy(self):
        edges = []
        for prot_id, features in uniprot_client.um.features.items():
            prot_node = ('UP', prot_id)
            for feature in features:
                feat_node = self.label('UPPRO', feature.id)
                edges.append((feat_node, prot_node,
                              {'type': 'partof'}))
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
            (self.label('INDRA_ACTIVITIES', source),
             self.label('INDRA_ACTIVITIES', target),
             {'type': 'isa'})
            for source, target in rels
        ]
        )

    def add_modification_hierarchy(self):
        self.add_edges_from([
            (self.label('INDRA_MODS', source),
             self.label('INDRA_MODS', 'modification'),
             {'type': 'isa'})
            for source in modtype_conditions
            if source != 'modification'
        ]
        )


def load_bio_ontology(reload=False):
    if reload or not os.path.exists(CACHE_FILE):
        ont = BioOntology()
        ont.initialize()
        # Try to create the folder first, if it fails, we don't cache
        if not os.path.exists(CACHE_DIR):
            try:
                os.makedirs(CACHE_DIR)
            except Exception:
                logger.warning('%s could not be created.' % CACHE_DIR)
                return ont
        # Try to dump the file next, if it fails, we don't cache
        try:
            logger.info('Caching INDRA bio ontology at %s' % CACHE_FILE)
            with open(CACHE_FILE, 'wb') as fh:
                pickle.dump(ont, fh, pickle.HIGHEST_PROTOCOL)
        except Exception:
            logger.warning('Failed to cache ontology at %s.' % CACHE_FILE)
        return ont
    else:
        logger.info('Loading INDRA bio ontology from cache at %s' % CACHE_FILE)
        with open(CACHE_FILE, 'rb') as fh:
            ont = pickle.load(fh)
            return ont


bio_ontology = load_bio_ontology(reload=False)
