from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from indra.statements import Phosphorylation, Agent, Evidence

logger = logging.getLogger('sparser')


class SparserProcessor(object):
    def __init__(self, xml_etree):
        self.tree = xml_etree
        self.statements = []

    def get_phosphorylations(self):
        phos_events = self.tree.findall("sem/ref/[@category='phosphorylate']")
        for event in phos_events:
            enzyme = event.find("var/[@name='agent']/ref")
            if enzyme is None:
                enz = None
            else:
                enz_name = enzyme.attrib.get('name')
                enz = Agent(enz_name)
            substrate = event.find("var/[@name='substrate']/ref")
            sub_name = substrate.attrib.get('name')
            sub = Agent(sub_name)
            st = Phosphorylation(enz, sub)
            self.statements.append(st)

    def _get_evidence(self, sem):
        ev = Evidence(source_api='sparser')
        return ev
