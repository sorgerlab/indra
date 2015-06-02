import sys
from pysb import *
from pysb.core import SelfExporter, InvalidComponentNameError, \
                      ComplexPattern, ReactionPattern

from belpy.statements import *

import xml.etree.ElementTree as ET

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
    def get_phosphorylation(self):
        phosphorylation_events = self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        for event in phosphorylation_events:
            sentence = self.get_text(event)
            agent = event.find(".//*[@role=':AGENT']")
            agent_name = self.get_text(agent)
            affected = event.find(".//*[@role=':AFFECTED']")
            affected_name = self.get_text(affected)
            mod = 'Phosphorylation'
            # TODO: query site one the TRIPS output is fixed
            site = ""
            self.belpy_stmts.append(Phosphorylation(agent_name,affected_name,mod,site,sentence,"",sentence,""))
            print agent_name, affected_name
        
       
    def make_model(self):
        model = Model()
        for stmt in self.belpy_stmts:
            stmt.monomers(model)
        for stmt in self.belpy_stmts:
            stmt.assemble(model)
        return model

if __name__ == '__main__':
    input_fname = 'phosphorylate.xml'
    if len(sys.argv)>1:
        input_fname = sys.argv[1]
    try:
        fh = open(input_fname,'rt')
    except IOError:
        print 'Could not open file %s' % input_fname
        exit()
    xml_string = fh.read()
    tp = TripsProcessor(xml_string)
    tp.get_complexes()
    tp.get_phosphorylation()
    model = tp.make_model()
