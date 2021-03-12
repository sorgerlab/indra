import os
import pickle
import logging
from indra.config import get_config
from ..ontology_graph import IndraOntology
from indra.util import read_unicode_csv
from indra.statements import modtype_conditions
from indra.resources import get_resource_path
from indra.statements.validate import assert_valid_db_refs


logger = logging.getLogger(__name__)


class BioOntology(IndraOntology):
    """Represents the ontology used for biology applications."""
    # The version is used to determine if the cached pickle is still valid
    # or not. When updating relevant resource files in INDRA, this version
    # should be incremented to "force" rebuilding the ontology to be consistent
    # with the underlying resource files.
    name = 'bio'
    version = '1.9'

    def __init__(self):
        super().__init__()

    def initialize(self, rebuild=False):
        if rebuild or not os.path.exists(CACHE_FILE):
            logger.info('Initializing INDRA bio ontology for the first time, '
                        'this may take a few minutes...')
            self._build()
            # Try to create the folder first, if it fails, we don't cache
            if not os.path.exists(CACHE_DIR):
                try:
                    os.makedirs(CACHE_DIR)
                except Exception:
                    logger.warning('%s could not be created.' % CACHE_DIR)
            # Try to dump the file next, if it fails, we don't cache
            try:
                logger.info('Caching INDRA bio ontology at %s' % CACHE_FILE)
                with open(CACHE_FILE, 'wb') as fh:
                    pickle.dump(self, fh, pickle.HIGHEST_PROTOCOL)
            except Exception:
                logger.warning('Failed to cache ontology at %s.' % CACHE_FILE)
        else:
            logger.info(
                'Loading INDRA bio ontology from cache at %s' % CACHE_FILE)
            with open(CACHE_FILE, 'rb') as fh:
                self.__dict__.update(pickle.load(fh).__dict__)

    def _build(self):
        # Add all nodes with annotations
        logger.info('Adding nodes...')
        self.add_hgnc_nodes()
        self.add_uniprot_nodes()
        self.add_famplex_nodes()
        self.add_obo_nodes()
        self.add_mesh_nodes()
        self.add_ncit_nodes()
        self.add_uppro_nodes()
        self.add_mirbase_nodes()
        self.add_chembl_nodes()
        self.add_hms_lincs_nodes()
        self.add_drugbank_nodes()
        # Add xrefs
        logger.info('Adding xrefs...')
        self.add_hgnc_uniprot_xrefs()
        self.add_hgnc_entrez_xrefs()
        self.add_famplex_xrefs()
        self.add_chemical_xrefs()
        self.add_ncit_xrefs()
        self.add_mesh_xrefs()
        self.add_mirbase_xrefs()
        self.add_hms_lincs_xrefs()
        self.add_biomappings()
        # Add hierarchies
        logger.info('Adding hierarchy...')
        self.add_famplex_hierarchy()
        self.add_obo_hierarchies()
        self.add_mesh_hierarchy()
        self.add_activity_hierarchy()
        self.add_modification_hierarchy()
        self.add_uppro_hierarchy()
        self._initialized = True
        # Build name to ID lookup
        logger.info('Building name lookup...')
        self._build_name_lookup()
        logger.info('Finished initializing bio ontology...')

    def add_hgnc_nodes(self):
        from indra.databases import hgnc_client
        withdrawns = set(hgnc_client.hgnc_withdrawn)
        nodes = [(self.label('HGNC', hid),
                  {'name': hname, 'obsolete': (hid in withdrawns)})
                 for (hid, hname) in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        from indra.databases import uniprot_client
        nodes = [(self.label('UP', uid), {'name': uname})
                 for (uid, uname)
                 in uniprot_client.um.uniprot_gene_name.items()]
        self.add_nodes_from(nodes)

    def add_uppro_nodes(self):
        from indra.databases import uniprot_client
        nodes = []
        for prot_id, features in uniprot_client.um.features.items():
            for feature in features:
                if feature.id is None:
                    continue
                node = self.label('UPPRO', feature.id)
                data = {'name': feature.name}
                nodes.append((node, data))
        self.add_nodes_from(nodes)

    def add_hgnc_uniprot_xrefs(self):
        from indra.databases import hgnc_client
        from indra.databases import uniprot_client
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

    def add_hgnc_entrez_xrefs(self):
        from indra.databases import hgnc_client
        edges = []
        for hid, eid in hgnc_client.entrez_ids.items():
            edges.append((self.label('HGNC', hid), self.label('EGID', eid),
                          {'type': 'xref', 'source': 'hgnc'}))
            edges.append((self.label('EGID', eid), self.label('HGNC', hid),
                          {'type': 'xref', 'source': 'hgnc'}))
        self.add_edges_from(edges)

    def add_famplex_nodes(self):
        nodes = []
        for row in read_unicode_csv(get_resource_path(
                os.path.join('famplex', 'entities.csv')), delimiter=','):
            entity = row[0]
            nodes.append((self.label('FPLX', entity),
                          {'name': entity}))
        self.add_nodes_from(nodes)

    def add_famplex_hierarchy(self):
        from indra.databases import hgnc_client
        edges = []
        for row in read_unicode_csv(get_resource_path(
                os.path.join('famplex', 'relations.csv')), delimiter=','):
            ns1, id1, rel, ns2, id2 = row
            if ns1 == 'HGNC':
                id1 = hgnc_client.get_hgnc_id(id1)
            edges.append((self.label(ns1, id1),
                          self.label(ns2, id2),
                          {'type': rel}))
        self.add_edges_from(edges)

    def add_famplex_xrefs(self):
        edges = []
        include_refs = {'PF', 'IP', 'GO', 'NCIT', 'ECCODE', 'HGNC_GROUP',
                        'MESH'}
        for row in read_unicode_csv(get_resource_path('famplex_map.tsv'),
                                    delimiter='\t'):
            ref_ns, ref_id, fplx_id = row
            if ref_ns not in include_refs:
                continue
            edges.append((self.label(ref_ns, ref_id),
                          self.label('FPLX', fplx_id),
                          {'type': 'xref', 'source': 'fplx'}))
            edges.append((self.label('FPLX', fplx_id),
                          self.label(ref_ns, ref_id),
                          {'type': 'xref', 'source': 'fplx'}))
        self.add_edges_from(edges)

    def add_obo_nodes(self):
        from indra.databases import obo_client
        namespaces = ['go', 'efo', 'hp', 'doid', 'chebi']
        nodes = []
        for ns in namespaces:
            oc = obo_client.OboClient(prefix=ns)
            for db_id, entry in oc.entries.items():
                label = self.label(ns.upper(), db_id)
                # There may be more general ways of doing this but this
                # handles the special case of BFO nodes in the EFO
                # ontology.
                if ns == 'efo' and db_id.startswith('BFO'):
                    label = db_id
                nodes.append((label,
                              {'name': entry['name']}))
        self.add_nodes_from(nodes)

    def add_obo_hierarchies(self):
        from indra.databases import obo_client
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
                    source_label = self.label(ns.upper(), db_id)
                    if ns == 'efo' and db_id.startswith('BFO'):
                        source_label = db_id
                    for target in targets:
                        target_label = self.label(ns.upper(), target)
                        if ns == 'efo' and target.startswith('BFO'):
                            target_label = target
                        edges.append((source_label,
                                      target_label,
                                      {'type': mapped_rel}))
        self.add_edges_from(edges)

    def add_chemical_xrefs(self):
        from indra.databases import chebi_client, drugbank_client
        mappings = [
            (chebi_client.chebi_chembl, 'CHEBI', 'CHEMBL', True),
            (chebi_client.chebi_pubchem, 'CHEBI', 'PUBCHEM', False),
            (chebi_client.pubchem_chebi, 'PUBCHEM', 'CHEBI', False),
            (chebi_client.hmdb_chebi, 'HMDB', 'CHEBI', True),
            (chebi_client.cas_chebi, 'CAS', 'CHEBI', True),
            (drugbank_client.drugbank_to_db, 'DRUGBANK', None, False),
            (drugbank_client.db_to_drugbank, None, 'DRUGBANK', False),

        ]
        edges = []
        data = {'type': 'xref', 'source': 'chebi'}

        def label_fix(ns, id):
            if ns == 'CHEBI' and not id.startswith('CHEBI'):
                id = 'CHEBI:%s' % id
            return self.label(ns, id)

        for map_dict, from_ns, to_ns, symmetric in mappings:
            for from_id, to_id in map_dict.items():
                # Here we assume if no namespace is given, then
                # we're dealing with a (namespace, id) tuple
                if from_ns is None:
                    from_ns_, from_id = from_id
                    to_ns_ = to_ns
                elif to_ns is None:
                    from_id, to_ns_ = from_id
                    from_ns_ = from_ns
                else:
                    from_ns_, to_ns_ = from_ns, to_ns
                source = label_fix(from_ns_, from_id)
                target = label_fix(to_ns_, to_id)
                edges.append((source, target, data))
                if symmetric:
                    edges.append((target, source, data))
        self.add_edges_from(edges)

    def add_mesh_nodes(self):
        from indra.databases import mesh_client
        nodes = [(self.label('MESH', mesh_id),
                  {'name': name})
                 for mesh_id, name in
                 mesh_client.mesh_id_to_name.items()]
        self.add_nodes_from(nodes)

    def add_mesh_xrefs(self):
        from indra.databases import mesh_client
        edges = []
        data = {'type': 'xref', 'source': 'gilda'}
        for mesh_id, (db_ns, db_id) in mesh_client.mesh_to_db.items():
            edges.append((self.label('MESH', mesh_id),
                          self.label(db_ns, db_id),
                          data))
        for (db_ns, db_id), mesh_id in mesh_client.db_to_mesh.items():
            edges.append((self.label(db_ns, db_id),
                          self.label('MESH', mesh_id),
                          data))
        self.add_edges_from(edges)

    def add_mesh_hierarchy(self):
        from indra.databases import mesh_client
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

    def add_biomappings(self):
        biomappings_tsv = get_resource_path('biomappings.tsv')
        edges = []
        for source_ns, source_id, _, target_ns, target_id, _ in \
                read_unicode_csv(biomappings_tsv, delimiter='\t'):
            edges.append((self.label(source_ns, source_id),
                          self.label(target_ns, target_id),
                          {'type': 'xref', 'source': 'biomappings'}))
            edges.append((self.label(target_ns, target_id),
                          self.label(source_ns, source_id),
                          {'type': 'xref', 'source': 'biomappings'}))
        self.add_edges_from(edges)

    def add_ncit_nodes(self):
        from indra.sources.trips.processor import ncit_map
        nodes = [(self.label('NCIT', ncit_id), {}) for ncit_id in ncit_map]
        self.add_nodes_from(nodes)

    def add_ncit_xrefs(self):
        from indra.sources.trips.processor import ncit_map
        edges = []
        for ncit_id, (target_ns, target_id) in ncit_map.items():
            edges.append((self.label('NCIT', ncit_id),
                          self.label(target_ns, target_id),
                          {'type': 'xref', 'source': 'ncit'}))
        self.add_edges_from(edges)

    def add_uppro_hierarchy(self):
        from indra.databases import uniprot_client
        edges = []
        for prot_id, features in uniprot_client.um.features.items():
            prot_node = self.label('UP', prot_id)
            for feature in features:
                if feature.id is None:
                    continue
                feat_node = self.label('UPPRO', feature.id)
                edges.append((feat_node, prot_node,
                              {'type': 'partof'}))
        self.add_edges_from(edges)

    def add_mirbase_nodes(self):
        from indra.databases import mirbase_client
        nodes = []
        for mirbase_id, name in mirbase_client._mirbase_id_to_name.items():
            nodes.append((self.label('MIRBASE', mirbase_id),
                          {'name': name}))
        self.add_nodes_from(nodes)

    def add_mirbase_xrefs(self):
        from indra.databases import mirbase_client
        edges = []
        for mirbase_id, hgnc_id in \
                mirbase_client._mirbase_id_to_hgnc_id.items():
            edges.append((self.label('MIRBASE', mirbase_id),
                          self.label('HGNC', hgnc_id),
                          {'type': 'xref', 'source': 'mirbase'}))
        for hgnc_id, mirbase_id in \
                mirbase_client._hgnc_id_to_mirbase_id.items():
            edges.append((self.label('HGNC', hgnc_id),
                          self.label('MIRBASE', mirbase_id),
                          {'type': 'xref', 'source': 'mirbase'}))
        self.add_edges_from(edges)

    def add_chembl_nodes(self):
        from indra.databases import chembl_client
        nodes = []
        for chembl_id, chembl_name in chembl_client.chembl_names.items():
            nodes.append((self.label('CHEMBL', chembl_id),
                          {'name': chembl_name}))
        self.add_nodes_from(nodes)

    def add_hms_lincs_nodes(self):
        from indra.databases.lincs_client import LincsClient
        lc = LincsClient()

        nodes = []
        for hmsl_id, data in lc._sm_data.items():
            if '-' in hmsl_id:
                hmsl_base_id, suffix = hmsl_id.split('-')
            else:
                hmsl_base_id, suffix = hmsl_id, None
            if suffix == '999':
                continue
            nodes.append((self.label('HMS-LINCS', hmsl_base_id),
                          {'name': data['Name']}))
        self.add_nodes_from(nodes)

    def add_hms_lincs_xrefs(self):
        from indra.databases.lincs_client import LincsClient
        lc = LincsClient()

        edges = []
        for hmsl_id, data in lc._sm_data.items():
            if '-' in hmsl_id:
                hmsl_base_id, suffix = hmsl_id.split('-')
            else:
                hmsl_base_id, suffix = hmsl_id, None
            if suffix == '999':
                continue
            refs = lc.get_small_molecule_refs(hmsl_id)
            for ref_ns, ref_id in refs.items():
                edges.append((self.label('HMS-LINCS', hmsl_base_id),
                              self.label(ref_ns, ref_id),
                              {'type': 'xref', 'source': 'hms-lincs'}))
                edges.append((self.label(ref_ns, ref_id),
                              self.label('HMS-LINCS', hmsl_base_id),
                              {'type': 'xref', 'source': 'hms-lincs'}))
        self.add_edges_from(edges)

    def add_drugbank_nodes(self):
        from indra.databases import drugbank_client
        nodes = []
        for db_id, db_name in drugbank_client.drugbank_names.items():
            nodes.append((self.label('DRUGBANK', db_id),
                          {'name': db_name}))
        self.add_nodes_from(nodes)

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

    def add_nodes_from(self, nodes_for_adding, **attr):
        for label, _ in nodes_for_adding:
            self.assert_valid_node(label)
        super().add_nodes_from(nodes_for_adding, **attr)

    def add_edges_from(self, ebunch_to_add, **attr):
        for source, target, _ in ebunch_to_add:
            for label in [source, target]:
                self.assert_valid_node(label)
        super().add_edges_from(ebunch_to_add, **attr)

    def assert_valid_node(self, label):
        db_ns, db_id = self.get_ns_id(label)
        if db_ns in {'INDRA_ACTIVITIES', 'INDRA_MODS'}:
            return
        try:
            assert_valid_db_refs({db_ns: db_id})
        except Exception as e:
            logger.warning('Invalid node: %s' % e)


CACHE_DIR = os.path.join((get_config('INDRA_RESOURCES') or
                          os.path.join(os.path.expanduser('~'), '.indra')),
                         '%s_ontology' % BioOntology.name,
                         BioOntology.version)
CACHE_FILE = os.path.join(CACHE_DIR, 'bio_ontology.pkl')
