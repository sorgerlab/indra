from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import collections
from indra.util import read_unicode_csv
from indra.literature import id_lookup
from indra.statements import *
from indra.databases import uniprot_client, hgnc_client

logger = logging.getLogger('sparser')


class SparserProcessor(object):
    def __init__(self, xml_etree):
        self.tree = xml_etree
        self.statements = []
        # Extract all sems by category
        self._sems = collections.defaultdict(list)
        for interp in self.tree.findall('interpretation'):
            sentence = interp.find('sentence-text').text
            sems = interp.findall('sem')
            for sem in sems:
                ref = sem.find('ref')
                if ref is not None:
                    category = ref.attrib['category']
                    self._sems[category].append((sem, sentence))
        # Get citation info
        pmcid = self.tree.attrib.get('pmcid')
        pmid = self.tree.attrib.get('pmid')
        self.pmid = None
        if pmid:
            self.pmid = pmid
        elif pmcid:
            ids = id_lookup(pmcid, 'pmcid')
            pmid = ids.get('pmid')
            if pmid is not None:
                self.pmid = pmid

    def get_modifications(self):
        mod_events = {}
        for mod_type in _mod_class_map.keys():
            mod_events[mod_type] = self._sems.get(mod_type)
        if not mod_events:
            return
        for mod_type, sems in mod_events.items():
            if not sems:
                continue
            indra_class = _mod_class_map.get(mod_type)
            if indra_class is None:
                logger.warning('Unhandled modification type %s' % mod_type)
                continue
            for event, sentence in sems:
                # Get enzyme agent
                enzyme = event.find("ref/var/[@name='agent']/ref")
                if enzyme is None:
                    enz = None
                else:
                    enz = self._get_agent_from_ref(enzyme)

                # Get substrate agent
                substrate = event.find("ref/var/[@name='substrate']/ref")
                # TODO: handle agent-or-substrate
                if substrate is None:
                    logger.debug('Skipping phosphorylation without substrate.')
                    continue
                sub = self._get_agent_from_ref(substrate)
                if sub is None:
                    logger.debug('Skipping phosphorylation without substrate.')
                    continue

                # Get site
                residue = None
                position = None
                site = event.find("ref/var/[@name='site']")
                if site is not None:
                    residue, position = self._get_site(site)

                # Get evidence
                ev = self._get_evidence(sentence)

                st = indra_class(enz, sub, residue, position, evidence=[ev])
                self.statements.append(st)

    def get_activations(self):
        act_events = self._sems.get('bio-activate')
        inh_events = self._sems.get('inhibit')
        for event_type, events in zip(['act', 'inh'],
                                      [act_events, inh_events]):
            if not events:
                continue
            for event, sentence in events:
                # Get the subject of the activation
                subj_ref = event.find("ref/var/[@name='by-means-of-or-agent']/ref")
                if subj_ref is None:
                    subj_ref = event.find("ref/var/[@name='agent']/ref")
                if subj_ref is None:
                    logger.debug('Skipping activation without subject.')
                    continue
                subj = self._get_agent_from_ref(subj_ref)
                if subj is None:
                    logger.debug('Skipping activation without subject.')
                    continue
                # Get the object of the activation
                obj_ref = event.find("ref/var/[@name='object']/ref")
                if obj_ref is None:
                    logger.debug('Skipping activation without object.')
                    continue
                obj = self._get_agent_from_ref(obj_ref)
                if obj is None:
                    logger.debug('Skipping activation without object.')
                    continue
                # Get evidence
                ev = self._get_evidence(sentence)
                if event_type == 'act':
                    st = Activation(subj, obj, evidence=[ev])
                else:
                    st = Inhibition(subj, obj, evidence=[ev])
                self.statements.append(st)

    def _get_agent_from_ref(self, ref):
        # TODO: handle collections
        if ref.attrib.get('category') == 'collection':
            #logger.warning('Skipping collection Agent.')
            return None

        # Find the name, uid and raw-text tags first and get their text
        # content if available
        uid_tag = ref.find("var/[@name='uid']")
        name_tag = ref.find("var/[@name='name']")
        text_tag = ref.find("var/[@name='raw-text']")
        if name_tag is not None and name_tag.text:
            name = name_tag.text
        else:
            name = None
        if uid_tag is not None and uid_tag.text:
            uid = uid_tag.text
        else:
            uid = None
        if text_tag is not None and text_tag.text:
            raw_text = text_tag.text
        else:
            raw_text = None

        db_refs = {}
        # Save raw text if available
        if raw_text:
            db_refs['TEXT'] = raw_text
        agent_name = raw_text
        # If we have a proper UID then we try to reconstruct an Agent from that
        if uid is not None and len(uid.split(':')) == 2:
            db_ns, db_id = uid.split(':')
            be_id = bioentities_map.get((db_ns, db_id))
            if be_id:
                db_refs[db_ns] = db_id
                db_refs['BE'] = be_id
                agent_name = be_id
            elif db_ns == 'UP':
                db_refs['UP'] = db_id
                gene_name = uniprot_client.get_gene_name(db_id)
                if gene_name:
                    agent_name = gene_name
                    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                    if hgnc_id:
                        db_refs['HGNC'] = hgnc_id
            elif db_ns == 'NCIT':
                db_refs['NCIT'] = db_id
                target = ncit_map.get(db_id)
                if target:
                    db_refs[target[0]] = target[1]
                    if target[0] == 'HGNC':
                        up_id = hgnc_client.get_uniprot_id(target[1])
                        agent_name = hgnc_client.get_hgnc_name(target[1])
                        if up_id:
                            db_refs['UP'] = up_id
                    elif target[0] == 'UP':
                        agent_name = uniprot_client.get_gene_name(target[1])
                        if agent_name:
                            hgnc_id = hgnc_client.get_hgnc_id(agent_name)
                            if hgnc_id:
                                db_refs['HGNC'] = hgnc_id
            elif db_ns in ['GO', 'CHEBI', 'FA', 'XFAM', 'MESH']:
                db_refs[db_ns] = db_id
            else:
                logger.warning('Unknown database name space %s' % db_ns)
        if not agent_name:
            if raw_text is not None:
                agent_name = raw_text
            else:
                return None

        assert(agent_name)

        agent = Agent(agent_name, db_refs=db_refs)
        return agent

    @staticmethod
    def _get_amino_acid(ref):
        name_tag = ref.find("var/[@name='name']")
        if name_tag is not None:
            return name_tag.text
        return None

    @staticmethod
    def _get_site(site):
        residue = None
        position = None
        residue_tag = site.find("ref/[@category='residue-on-protein']")
        if residue_tag is not None:
            residue_aa_tag = residue_tag.find("var/[@name='amino-acid']/ref")
            if residue_aa_tag is not None:
                residue = SparserProcessor._get_amino_acid(residue_aa_tag)
                position_tag = residue_tag.find("var/[@name='position']")
                if position_tag is not None:
                    position = position_tag.text.strip()
        if residue and residue.startswith('phospho'):
            residue = residue[6:]
            try:
                residue = get_valid_residue(residue)
            except InvalidResidueError:
                logger.error('Invalid residue found: %s' % residue)
                residue = None
        return residue, position

    def _get_evidence(self, text):
        ev = Evidence(source_api='sparser', pmid=self.pmid, text=text)
        return ev

_mod_class_map = {
    'phosphorylate': Phosphorylation,
    'dephosphorylate': Dephosphorylation,
    'ubiquitination': Ubiquitination,
    'deubiquitination': Deubiquitination,
    'farnesylation': Farnesylation,
    'defarnesylation': Defarnesylation
    }

def _read_ncit_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../resources/ncit_map.tsv')
    ncit_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    next(csv_rows)
    for row in csv_rows:
        ncit_id = row[0]
        target_ns = row[1]
        target_id = row[2]
        ncit_map[ncit_id] = (target_ns, target_id)
    return ncit_map

ncit_map = _read_ncit_map()

def _read_bioentities_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../resources/bioentities_map.tsv')
    bioentities_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        source_ns = row[0]
        source_id = row[1]
        be_id = row[2]
        bioentities_map[(source_ns, source_id)] = be_id
    return bioentities_map

bioentities_map = _read_bioentities_map()
