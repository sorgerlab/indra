import re
import warnings
import pickle

import xml.etree.ElementTree as ET

from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client


residue_names = {
    'S': 'Serine',
    'T': 'Threonine',
    'Y': 'Tyrosine',
    'SER': 'Serine',
    'THR': 'Threonine',
    'TYR': 'Tyrosine',
    'SERINE': 'Serine',
    'THREONINE': 'Threonine',
    'TYROSINE': 'Tyrosine'
    }


mod_names = {
    'PHOSPHORYLATION': 'Phosphorylation'
    }


class TripsProcessor(object):
    def __init__(self, xml_string):
        self.tree = ET.fromstring(xml_string)
        self.statements = []
        self._hgnc_cache = self._load_hgnc_cache()
        self._static_events = self._find_static_events()

    def get_activations(self):
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        for event in act_events:
            sentence = self._get_text(event)

            # Get the activating agent in the event
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                continue
            agent_id = agent.attrib['id']
            agent_name = self._get_name_by_id(agent_id)
            agent_agent = Agent(agent_name)

            # Get the activated agent in the event
            affected = event.find(".//*[@role=':AFFECTED']")
            affected_id = affected.attrib['id']
            affected_name = self._get_name_by_id(affected_id)
            affected_agent = Agent(affected_name)

            ev = Evidence(source_api='trips', text=sentence)
            self.statements.append(ActivityActivity(agent_agent, 'act',
                                    'increases', affected_agent, 'act',
                                    evidence=ev))

    def get_activating_mods(self):
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        for event in act_events:
            if event.attrib['id'] in self._static_events:
                continue
            sentence = self._get_text(event)
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                msg = 'Skipping activation event with no affected term.'
                warnings.warn(msg)
                continue

            affected_id = affected.attrib['id']
            affected_name = self._get_name_by_id(affected_id)
            affected_agent = Agent(affected_name)
            precond_event_ref = \
                self.tree.find("TERM/[@id='%s']/features/inevent" % affected_id)
            if precond_event_ref is None:
                # This means that it is not an activating modification
                continue
            precond_id = precond_event_ref.find('event').attrib['id']
            precond_event = self.tree.find("EVENT[@id='%s']" % precond_id)
            mod, mod_pos = self._get_mod_site(precond_event)

            ev = Evidence(source_api='trips', text=sentence)
            self.statements.append(ActivityModification(affected_agent, mod,
                                    mod_pos, 'increases', 'Active',
                                    evidence=ev))

    def get_complexes(self):
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        for event in bind_events:
            if event.attrib['id'] in self._static_events:
                continue

            sentence = self._get_text(event)

            arg1 = event.find("arg1")
            if arg1 is None:
                msg = 'Skipping complex missing arg1.'
                warnings.warn(msg)
                continue
            agent1 = self._get_agent_by_id(arg1.attrib['id'], event.attrib['id'])

            arg2 = event.find("arg2")
            if arg2 is None:
                msg = 'Skipping complex missing arg2.'
                warnings.warn(msg)
                continue
            agent2 = self._get_agent_by_id(arg2.attrib['id'], event.attrib['id'])

            # Information on binding site is either attached to the agent term
            # in a features/site tag or attached to the event itself in 
            # a site tag
            site_feature = self._find_in_term(arg1.attrib['id'], 'features/site')
            if site_feature is not None:
                sites, positions = self._get_site_by_id(site_id)
                print sites, positions

            site_feature = self._find_in_term(arg2.attrib['id'], 'features/site')
            if site_feature is not None:
                sites, positions = self._get_site_by_id(site_id)
                print sites, positions
            
            site = event.find("site")
            if site is not None:
                sites, positions = self._get_site_by_id(site.attrib['id'])
                print sites, positions
            
            if agent1 is None or agent2 is None:
                warnings.warn('Complex with missing members')
                continue
            self.statements.append(Complex([agent1, agent2]))

    def get_phosphorylation(self):
        phosphorylation_events = \
            self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        for event in phosphorylation_events:
            if event.attrib['id'] in self._static_events:
                continue

            sentence = self._get_text(event)
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                warnings.warn('Skipping phosphorylation event with no agent.')
                continue
            if agent.find("type").text == 'ONT::MACROMOLECULAR-COMPLEX':
                complex_id = agent.attrib['id']
                complex_term = self.tree.find("TERM/[@id='%s']" % complex_id)
                components = complex_term.find("components")
                terms = components.findall("term")
                term_names = []
                for t in terms:
                    term_names.append(self._get_name_by_id(t.attrib['id']))
                agent_name = term_names[0]
                agent_bound = Agent(term_names[1])
                agent_agent = Agent(agent_name,
                            bound_conditions=[BoundCondition(agent_bound,True)])
            else:
                agent_agent = self._get_agent_by_id(agent.attrib['id'],
                                                    event.attrib['id'])
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                warnings.warn('Skipping phosphorylation event with no '
                              'affected term.')
                continue
            affected_agent = self._get_agent_by_id(affected.attrib['id'],
                                                   event.attrib['id'])
            mod, mod_pos = self._get_mod_site(event)
            # TODO: extract more information about text to use as evidence
            ev = Evidence(source_api='trips', text=sentence)
            # Assuming that multiple modifications can only happen in
            # distinct steps, we add a statement for each modification
            # independently

            mod_types = event.findall('predicate/mods/mod/type')
            # Transphosphorylation
            if 'ONT::ACROSS' in [mt.text for mt in mod_types]:
                agent_bound = Agent(affected_agent.name)
                agent_agent.bound_conditions = \
                                           [BoundCondition(agent_bound, True)]
                for m, p in zip(mod, mod_pos):
                    self.statements.append(Transphosphorylation(agent_agent,
                                        m, p, evidence=ev))
            # Dephosphorylation
            elif 'ONT::MANNER-UNDO' in [mt.text for mt in mod_types]:
                for m, p in zip(mod, mod_pos):
                    self.statements.append(Dephosphorylation(agent_agent,
                                        affected_agent, m, p, evidence=ev))
            # Autophosphorylation
            elif agent.attrib['id'] == affected.attrib['id']:
                for m, p in zip(mod, mod_pos):
                    self.statements.append(Autophosphorylation(agent_agent,
                                        m, p, evidence=ev))
            # Regular phosphorylation
            else:
                for m, p in zip(mod, mod_pos):
                    self.statements.append(Phosphorylation(agent_agent,
                                            affected_agent, m, p, evidence=ev))

    def _get_agent_by_id(self, entity_id, event_id):
        term = self.tree.find("TERM/[@id='%s']" % entity_id)
        if term is None:
            return None

        # Extract database references
        try:
            dbid = term.attrib["dbid"]
            dbids = dbid.split('|')
            db_refs_dict = dict([d.split(':') for d in dbids])
        except KeyError:
            db_refs_dict = {}

        # If the entity is a complex
        if term.find("type").text == 'ONT::MACROMOLECULAR-COMPLEX':
            complex_id = entity_id
            complex_term = self.tree.find("TERM/[@id='%s']" % complex_id)
            components = complex_term.find("components")
            if components is None:
                warnings.warn('Complex without components')
                return None
            terms = components.findall('component')
            term_names = []
            agents = []
            for t in terms:
                agents.append(self._get_agent_by_id(t.attrib['id'], None))
            # We assume that the first agent mentioned in the description of
            # the complex is the one that mediates binding
            agent = agents[0]
            agent.bound_conditions = \
                            [BoundCondition(ag, True) for ag in agents[1:]]
        # If the entity is not a complex
        else:
            agent_name = self._get_name_by_id(entity_id)
            agent = Agent(agent_name, db_refs=db_refs_dict)
            precond_event_ref = \
                self.tree.find("TERM/[@id='%s']/features/inevent" % entity_id)
            # Extract preconditions of the agent
            if precond_event_ref is not None:
                # Find the event describing the precondition
                precond_id = precond_event_ref.find('event').attrib['id']
                precond_event = self.tree.find("EVENT[@id='%s']" % precond_id)
                if precond_id == event_id:
                    warnings.warn('Circular reference to event %s.' %
                                  precond_id)
                else:
                    precond_event_type = precond_event.find('type').text
                    # Binding precondition
                    if precond_event_type == 'ONT::BIND':
                        arg1 = precond_event.find('arg1')
                        arg2 = precond_event.find('arg2')
                        mod = precond_event.findall('mods/mod')
                        if arg1 is None:
                            arg2_name = self._get_name_by_id(arg2.attrib['id'])
                            bound_agent = Agent(arg2_name)
                        elif arg2 is None:
                            arg1_name = self._get_name_by_id(arg1.attrib['id'])
                            bound_agent = Agent(arg1_name)
                        else:
                            arg1_name = self._get_name_by_id(arg1.attrib['id'])
                            arg2_name = self._get_name_by_id(arg2.attrib['id'])
                            if arg1_name == agent_name:
                                bound_agent = Agent(arg2_name)
                            else:
                                bound_agent = Agent(arg1_name)
                        # Look for negative flag either in precondition event
                        # predicate tag or in the term itself
                        # (after below, neg_flag will be an object, or None)
                        neg_flag = precond_event.find(
                                        'predicate/mods/mod[type="ONT::NEG"]')
                        negation_sign = precond_event.find('negation')
                        if negation_sign is not None:
                            if negation_sign.text == '+':
                                neg_flag = True
                        # (after this, neg_flag will be a boolean value)
                        neg_flag = neg_flag or \
                                   term.find('mods/mod[type="ONT::NEG"]')
                        negation_sign = precond_event.find('predicate/negation')
                        if negation_sign is not None:
                            if negation_sign.text == '+':
                                neg_flag = True

                        if neg_flag:
                            bc = BoundCondition(bound_agent, False)
                        else:
                            bc = BoundCondition(bound_agent, True)
                        agent.bound_conditions = [bc]

                    # Phosphorylation precondition
                    elif precond_event_type == 'ONT::PHOSPHORYLATION':
                        mod, mod_pos = self._get_mod_site(precond_event)
                        for m, mp in zip(mod, mod_pos):
                            agent.mods.append(m)
                            agent.mod_sites.append(mp)
            
        return agent

    def _find_in_term(self, term_id, path):
        tag = self.tree.find("TERM[@id='%s']/%s" % (term_id, path))
        return tag
    
    def _get_text(self, element):
        text_tag = element.find("text")
        text = text_tag.text
        return text

    def _get_hgnc_name(self, hgnc_id):
        try:
            hgnc_name = self._hgnc_cache[hgnc_id]
        except KeyError:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)

            self._hgnc_cache[hgnc_id] = hgnc_name
        return hgnc_name

    def _get_valid_component_name(self, name):
        name = name.replace('-', '_')
        return name

    def _get_name_by_id(self, entity_id):
        entity_term = self.tree.find("TERM/[@id='%s']" % entity_id)
        name = entity_term.find("name")
        if name is None:
            warnings.warn('Entity without a name')
            return ''
        try:
            dbid = entity_term.attrib["dbid"]
        except:
            warnings.warn('No grounding information for %s' % name.text)
            return self._get_valid_component_name(name.text)
        dbids = dbid.split('|')
        hgnc_ids = [i for i in dbids if i.startswith('HGNC')]
        up_ids = [i for i in dbids if i.startswith('UP')]
        #TODO: handle protein families like 14-3-3 with IDs like
        # XFAM:PF00244.15, FA:00007
        if hgnc_ids:
            if len(hgnc_ids) > 1:
                warnings.warn('%d HGNC IDs reported.' % len(hgnc_ids))
            hgnc_id = re.match(r'HGNC\:([0-9]*)', hgnc_ids[0]).groups()[0]
            hgnc_name = self._get_hgnc_name(hgnc_id)
            return self._get_valid_component_name(hgnc_name)
        elif up_ids:
            if len(hgnc_ids) > 1:
                warnings.warn('%d UniProt IDs reported.' % len(up_ids))
            up_id = re.match(r'UP\:([A-Z0-9]*)', up_ids[0]).groups()[0]
            up_rdf = up_client.query_protein(up_id)
            # First try to get HGNC name
            hgnc_name = up_client.get_hgnc_name(up_rdf)
            if hgnc_name is not None:
                return self._get_valid_component_name(hgnc_name)
            # Next, try to get the gene name
            gene_name = up_client.get_gene_name(up_rdf)
            if gene_name is not None:
                return self._get_valid_component_name(gene_name)
        # By default, return the text of the name tag
        name_txt = name.text.strip('|')
        return self._get_valid_component_name(name_txt)

    # Get all the sites recursively based on a term id.
    def _get_site_by_id(self, site_id):
        all_residues = []
        all_pos = []
        site_term = self.tree.find("TERM/[@id='%s']" % site_id)
        aggregate = site_term.find('aggregate')
        if aggregate is not None:
            for member in aggregate.getchildren():
                residue, pos = self._get_site_by_id(member.attrib['id'])
                all_residues.extend(residue)
                all_pos.extend(pos)
        else:
            site_type = site_term.find("type").text
            site_name = site_term.find("name").text
            if site_type == 'ONT::MOLECULAR-SITE':
                residue = site_term.find('features/site/name').text.upper()
                pos = site_term.find('features/site/pos').text.upper()
            elif site_type == 'ONT::RESIDUE':
                # Example name: TYROSINE-RESIDUE
                residue = site_name.split('-')[0]
                pos = None
            elif site_type == 'ONT::AMINO-ACID':
                residue = site_name
                pos = None
            return (residue, ), (pos, )
        return all_residues, all_pos

    def _get_mod_site(self, event):
        mod_type = event.find('type')
        mod_type_name = mod_names[mod_type.text.split('::')[1]]

        site_tag = event.find("site")
        if site_tag is None:
            return [mod_type_name], ['']
        site_id = site_tag.attrib['id']
        residues, mod_pos = self._get_site_by_id(site_id)
        mod = [mod_type_name+residue_names[r] for r in residues]
        return mod, mod_pos

    def _find_static_events(self):
        inevent_tags = self.tree.findall("TERM/features/inevent/event")
        static_events = []
        for ie in inevent_tags:
            static_events.append(ie.attrib['id'])
        return static_events

    def _load_hgnc_cache(self):
        try:
            fh = open('hgnc_cache.pkl', 'rb')
        except IOError:
            return {}
        return pickle.load(fh)

    def _dump_hgnc_cache(self):
        with open('hgnc_cache.pkl', 'wb') as fh:
            pickle.dump(self._hgnc_cache, fh)

if __name__ == '__main__':
    tp = TripsProcessor(open('wchen-v3.xml', 'rt').read())
    tp._find_static_events()
    tp.get_complexes()
    tp.get_phosphorylation()
    tp.get_activating_mods()
