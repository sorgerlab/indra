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
                category = sem.find('ref').attrib['category']
                self._sems[category].append((sem, sentence))

    def get_phosphorylations(self):
        phos_events = self._sems.get('phosphorylate')
        for event, sentence in phos_events:
            # Get enzyme agent
            enzyme = event.find("ref/var/[@name='agent']/ref")
            if enzyme is None:
                enz = None
            else:
                enz = self._get_agent_from_ref(enzyme)

            # Get substrate agent
            substrate = event.find("ref/var/[@name='substrate']/ref")
            sub = self._get_agent_from_ref(substrate)

            # Get site
            residue = None
            position = None
            site = event.find("ref/var/[@name='site']")
            if site is not None:
                residue, position = self._get_site(site)

            st = Phosphorylation(enz, sub, residue, position)
            self.statements.append(st)

    def _get_agent_from_ref(self, ref):
        name = ref.attrib.get('name')
        if name is None:
            name = ref.find("var/[@name='name']").text
            uid = ref.find("var/[@name='uid']").text
        else:
            uid = ref.attrib.get('uid')
        db_refs = {}
        if uid.startswith('UP:'):
            up_mnemonic = uid[3:]
            up_id = uniprot_client.get_id_from_mnemonic(up_mnemonic)
            name = uniprot_client.get_gene_name(up_id)
            db_refs['UP'] = up_id

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
            residue = SparserProcessor._get_amino_acid(residue_aa_tag)
            position_tag = residue_tag.find("var/[@name='position']")
            if position_tag is not None:
                position = position_tag.text.strip()
        return residue, position

    def _get_evidence(self, sem):
        ev = Evidence(source_api='sparser')
        return ev
