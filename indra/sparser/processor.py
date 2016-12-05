from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import collections
from indra.literature import id_lookup
from indra.statements import *
from indra.databases import uniprot_client

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

    def _get_agent_from_ref(self, ref):
        # TODO: handle collections
        if ref.attrib.get('category') == 'collection':
            logger.warning('Skipping collection Agent.')
            return None
        name_tag = ref.find("var/[@name='name']")
        if name_tag is not None:
            name = name_tag.text
        else:
            return None
        uid_tag = ref.find("var/[@name='uid']")
        if uid_tag is not None:
            uid = uid_tag.text
        else:
            uid = None

        db_refs = {}
        if uid is not None and uid.startswith('UP:'):
            up_mnemonic = uid[3:]
            up_id = uniprot_client.get_id_from_mnemonic(up_mnemonic)
            if up_id is not None:
                up_name = uniprot_client.get_gene_name(up_id)
                if up_name is not None:
                    name = up_name
                db_refs['UP'] = up_id

        assert name is not None

        agent = Agent(name, db_refs=db_refs)
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