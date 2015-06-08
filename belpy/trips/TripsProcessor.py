from pysb import *
from pysb.core import SelfExporter, InvalidComponentNameError, \
                      ComplexPattern, ReactionPattern

from belpy.statements import *

import xml.etree.ElementTree as ET

residue_names = {
    'SER':'Serine',
    'THR':'Threonine',
    'TYR':'Tyrosine'
    }

class TripsProcessor(object):
    def __init__(self, xml_string):
        self.tree = ET.fromstring(xml_string)
        self.belpy_stmts = []

    def get_text(self,element):
        text_tag = element.find("text")
        text = text_tag.text
        return text

    def get_complexes(self):
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        for event in bind_events:
            sentence = self.get_text(event)
            arg1 = event.find("arg1")
            arg1_name = self.get_text(arg1)
            arg2 = event.find("arg2")
            arg2_name = self.get_text(arg2)
            self.belpy_stmts.append(Complex((arg1_name,arg2_name)))

    def get_name_by_id(self,entity_id):
        entity_term = self.tree.find("TERM/[@id='%s']" % entity_id)
        name = entity_term.find("name")
        return name.text

    # Get all the sites recursively based on a term id. 
    def get_site_by_id(self,site_id):
        all_residues = []
        all_pos = []
        site_term = self.tree.find("TERM/[@id='%s']" % site_id)
        subterms = site_term.find("subterms")
        print subterms
        if subterms is not None:
            for s in subterms.getchildren():
                residue, pos = self.get_site_by_id(s.text)
                all_residues.extend(residue)
                all_pos.extend(pos)
        else:
            site_name = site_term.find("name")
            # Example name: SER-222
            residue, pos = site_name.text.split('-')
            return (residue, ), (pos, )
        return all_residues, all_pos
        

    def get_phosphorylation(self):
        phosphorylation_events = self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        for event in phosphorylation_events:
            sentence = self.get_text(event)
            agent = event.find(".//*[@role=':AGENT']")
            agent_name = self.get_name_by_id(agent.attrib['id'])
            affected = event.find(".//*[@role=':AFFECTED']")
            affected_name = self.get_name_by_id(affected.attrib['id'])
            site_tag = event.find("site")
            site_id = site_tag.attrib['id']
            residues, mod_pos = self.get_site_by_id(site_id)
            mod = 'Phosphorylation'
            mod = ['Phospho'+residue_names[r] for r in residues]
            # TODO: extract more information about text to use as evidence
            citation = ''
            evidence = sentence
            annotations = ''
            # Assuming that multiple modifications can only happen in 
            # distinct steps, we add a statement for each modification
            # independently
            for m,p in zip(mod, mod_pos):
                self.belpy_stmts.append(Phosphorylation(agent_name,affected_name,
                                    m,p,sentence,citation,evidence,annotations))
    
if __name__ == '__main__':
    tp = TripsProcessor(open('phosphosite2.xml','rt').read())
    tp.get_phosphorylation()


