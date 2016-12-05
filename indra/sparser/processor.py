from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import collections
from indra.statements import Phosphorylation, Agent, Evidence
from indra.databases import uniprot_client

logger = logging.getLogger('sparser')


class SparserProcessor(object):
    def __init__(self, xml_etree):
        self.tree = xml_etree
        self.statements = []
        self._sems = collections.defaultdict(list)
        for interp in self.tree.findall('interpretation'):
            sentence = interp.find('sentence-text').text
            sems = interp.findall('sem')
            for sem in sems:
                ref = sem.find('ref')
                if ref is not None:
                    category = ref.attrib['category']
                    self._sems[category].append((sem, sentence))

    def get_phosphorylations(self):
        phos_events = self._sems.get('phosphorylate')
        if not phos_events:
            return
        for event, sentence in phos_events:
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
                continue
            sub = self._get_agent_from_ref(substrate)
            if sub is None:
                continue

            # Get site
            residue = None
            position = None
            site = event.find("ref/var/[@name='site']")
            if site is not None:
                residue, position = self._get_site(site)

            st = Phosphorylation(enz, sub, residue, position)
            self.statements.append(st)

    def _get_agent_from_ref(self, ref):
        # TODO: handle collections
        if ref.attrib.get('category') == 'collection':
            logger.warning('Skipping collection Agent.')
            return None
        name = ref.attrib.get('name')
        uid = ref.attrib.get('uid')
        if name is None:
            name_tag = ref.find("var/[@name='name']")
            if name_tag is not None:
                name = name_tag.text
            else:
                return None
        if uid is None:
            uid_tag = ref.find("var/[@name='uid']")
            if uid_tag is not None:
                uid = uid_tag.text

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

    def _get_evidence(self, sem):
        ev = Evidence(source_api='sparser')
        return ev
