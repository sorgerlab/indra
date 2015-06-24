import re
import warnings
import ipdb

import xml.etree.ElementTree as ET
from BeautifulSoup import BeautifulSoup

from belpy.statements import *
import belpy.databases.hgnc_client as hgnc_client
import belpy.databases.uniprot_client as up_client

residue_names = {
    'SER': 'Serine',
    'THR': 'Threonine',
    'TYR': 'Tyrosine'
    }

mod_names = {
    'PHOSPHORYLATION': 'Phosphorylation'
    }


class TripsProcessor(object):
    def __init__(self, xml_string):
        self.tree = ET.fromstring(xml_string)
        self.belpy_stmts = []
        self.hgnc_cache = {}

    def get_text(self, element):
        text_tag = element.find("text")
        text = text_tag.text
        return text
    
    def get_hgnc_name(self, hgnc_id):
        try:
            hgnc_name = self.hgnc_cache[hgnc_id]
        except KeyError:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
            self.hgnc_cache[hgnc_id] = hgnc_name
        return hgnc_name

    def get_name_by_id(self, entity_id):
        entity_term = self.tree.find("TERM/[@id='%s']" % entity_id)
        name = entity_term.find("name")
        try:
            dbid = entity_term.attrib["dbid"]
        except:
            warnings.warn('No grounding information for %s' % name)
            return name.text
        dbids = dbid.split('|')
        hgnc_ids = [i for i in dbids if i.startswith('HGNC')]
        up_ids = [i for i in dbids if i.startswith('UP')]

        if hgnc_ids:
            if len(hgnc_ids) > 1:
                warnings.warn('%d HGNC IDs reported.' % len(hgnc_ids))
            hgnc_id = re.match(r'HGNC\:([0-9]*)', hgnc_ids[0]).groups()[0]
            hgnc_name = self.get_hgnc_name(hgnc_id)
            return hgnc_name
        else:
            if up_ids:
                if len(hgnc_ids) > 1:
                    warnings.warn('%d UniProt IDs reported.' % len(up_ids))
                up_id = re.match(r'UP\:([A-Z0-9]*)', up_ids[0]).groups()[0]
                up_rdf = up_client.query_protein(up_id)
                # First try to get HGNC name
                hgnc_name = up_client.get_hgnc_name(up_rdf)
                if hgnc_name is not None:
                    return hgnc_name
                # Next, try to get the gene name
                gene_name = up_client.get_gene_name(up_rdf)
                if gene_name is not None:
                    return gene_name
        # By default, return the text of the name tag
        return name.text

    # Get all the sites recursively based on a term id.
    def get_site_by_id(self, site_id):
        all_residues = []
        all_pos = []
        site_term = self.tree.find("TERM/[@id='%s']" % site_id)
        subterms = site_term.find("subterms")
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

    def get_activating_mods(self):
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        for event in act_events:
            sentence = self.get_text(event)
            affected = event.find(".//*[@role=':AFFECTED']")
            affected_id = affected.attrib['id']
            affected_name = self.get_name_by_id(affected_id)
            precond_event_ref = self.tree.find("TERM/[@id='%s']/features/inevent"
                                               % affected_id)
            if precond_event_ref is None:
                warnings.warn('Skipping activation event with \
                              no precondition event.')
            precond_id = precond_event_ref.find('eventID').text
            precond_event = self.tree.find("EVENT[@id='%s']" % precond_id)
            mod, mod_pos = self.get_mod_site(precond_event)
            citation = ''
            evidence = sentence
            annotations = ''
            self.belpy_stmts.append(ActivityModification(affected_name, mod,
                                    mod_pos, 'DirectlyIncreases', 'Active',
                                    sentence, citation, evidence, annotations))

    def get_complexes(self):
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        for event in bind_events:
            sentence = self.get_text(event)
            arg1 = event.find("arg1")
            arg1_name = self.get_name_by_id(arg1.attrib['id'])
            arg2 = event.find("arg2")
            arg2_name = self.get_name_by_id(arg2.attrib['id'])
            self.belpy_stmts.append(Complex((arg1_name, arg2_name)))

    def get_mod_site(self, event):
        mod_type = event.find('type')
        mod_type_name = mod_names[mod_type.text.split('::')[1]]

        site_tag = event.find("site")
        if site_tag is None:
            return [mod_type_name], ['']
        site_id = site_tag.attrib['id']
        residues, mod_pos = self.get_site_by_id(site_id)
        mod = [mod_type_name+residue_names[r] for r in residues]
        return mod, mod_pos

    def get_phosphorylation(self):
        phosphorylation_events = self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        ipdb.set_trace()
        for event in phosphorylation_events:
            sentence = self.get_text(event)
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                warnings.warn('Skipping phosphorylation event with no agent.')
                continue
            agent_name = self.get_name_by_id(agent.attrib['id'])
            affected = event.find(".//*[@role=':AFFECTED']")
            affected_name = self.get_name_by_id(affected.attrib['id'])
            mod, mod_pos = self.get_mod_site(event)
            # TODO: extract more information about text to use as evidence
            citation = ''
            evidence = sentence
            annotations = ''
            # Assuming that multiple modifications can only happen in
            # distinct steps, we add a statement for each modification
            # independently

            for m, p in zip(mod, mod_pos):
                self.belpy_stmts.append(Phosphorylation(agent_name, affected_name,
                                        m, p, sentence, citation, evidence, annotations))

if __name__ == '__main__':
    tp = TripsProcessor(open('wchen-v3.xml', 'rt').read())
    tp.get_complexes()
    tp.get_phosphorylation()
    tp.get_activating_mods()
