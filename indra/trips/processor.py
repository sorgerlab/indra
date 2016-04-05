import re
import warnings
import pickle
import operator

import xml.etree.ElementTree as ET

from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client


mod_names = {
    'PHOSPHORYLATION': 'phosphorylation'
    }


class TripsProcessor(object):
    def __init__(self, xml_string):
        try:
            self.tree = ET.fromstring(xml_string)
        except ET.ParseError:
            print 'Could not parse XML string'
            self.tree = None
            return
        # Get the document ID from the EKB tag. This is the PMC ID when available.
        self.doc_id = self.tree.attrib.get('id')
        # Store all paragraphs and store all sentences in a data structure
        paragraph_tags = self.tree.findall('input/paragraphs/paragraph')
        sentence_tags = self.tree.findall('input/sentences/sentence')
        self.paragraphs = {p.attrib['id']: p.text for p in paragraph_tags}
        self.sentences = {s.attrib['id']: s.text for s in sentence_tags}

        self.statements = []
        self._static_events = self._find_static_events()
        self.get_all_events()
        self.extracted_events = {k:[] for k in self.all_events.keys()}
        print 'All events by type'
        print '------------------'
        for k, v in self.all_events.iteritems():
            print k, len(v)
        print '------------------'

    def get_all_events(self):
        self.all_events = {}
        events = self.tree.findall('EVENT')
        for e in events:
            event_id = e.attrib['id']
            if event_id in self._static_events:
                continue
            event_type = e.find('type').text
            try:
                self.all_events[event_type].append(event_id)
            except KeyError:
                self.all_events[event_type] = [event_id]

    def get_evidence_text(self, event_tag):
        '''
        Pieces of text linked to an EVENT are fragments of a sentence. The
        EVENT refers to the paragraph ID and the "uttnum", which corresponds 
        to a sentence ID. Here we find and return the full sentence from which 
        the event was taken.
        '''
        par_id = event_tag.attrib.get('paragraph')
        uttnum = event_tag.attrib.get('uttnum')
        event_text = event_tag.find('text')
        if self.sentences is not None and uttnum is not None:
            sentence = self.sentences[uttnum]
        elif event_text is not None:
            sentence = event_text.text
        else:
            sentence = None
        return sentence

    def get_activations(self):
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        inact_events = self.tree.findall("EVENT/[type='ONT::DEACTIVATE']")
        inact_events += self.tree.findall("EVENT/[type='ONT::INHIBIT']")
        for event in (act_events + inact_events):
            sentence = self.get_evidence_text(event)
            ev = Evidence(source_api='trips', text=sentence,
                          pmid=self.doc_id)

            # Get the activating agent in the event
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                continue
            agent_id = agent.attrib['id']
            agent_name = self._get_name_by_id(agent_id)
            if agent_name is None:
                warnings.warn(
                    'Skipping activation with missing activator agent')
                continue
            activator_agent = Agent(agent_name)

            # Get the activated agent in the event
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                warnings.warn(
                    'Skipping activation with missing affected agent')
                continue
            affected_id = affected.attrib['id']
            affected_name = self._get_name_by_id(affected_id)
            if affected_name is None:
                warnings.warn(
                    'Skipping activation with missing affected agent')
                continue

            affected_agent = Agent(affected_name)

            if event.find('type').text == 'ONT::ACTIVATE':
                rel = 'increases'
                activator_act = 'Activity'
                self.extracted_events['ONT::ACTIVATE'].append(event.attrib['id'])
            elif event.find('type').text == 'ONT::INHIBIT':
                rel = 'decreases'
                activator_act = None
                self.extracted_events['ONT::INHIBIT'].append(event.attrib['id'])
            elif event.find('type').text == 'ONT::DEACTIVATE':
                rel = 'decreases'
                activator_act = 'Activity'
                self.extracted_events['ONT::DEACTIVATE'].append(event.attrib['id'])

            self.statements.append(ActivityActivity(activator_agent, activator_act,
                                    rel, affected_agent, 'Activity',
                                    evidence=ev))

    def get_activating_mods(self):
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        for event in act_events:
            if event.attrib['id'] in self._static_events:
                continue
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                msg = 'Skipping activation event with no affected term.'
                warnings.warn(msg)
                continue

            affected_id = affected.attrib['id']
            affected_name = self._get_name_by_id(affected_id)
            if affected_name is None:
                warnings.warn(
                    'Skipping activating modification with missing' +\
                    'affected agent')
                continue
            affected_agent = Agent(affected_name)
            precond_event_ref = \
                self.tree.find("TERM/[@id='%s']/features/inevent" % affected_id)
            if precond_event_ref is None:
                # This means that it is not an activating modification
                continue
            precond_id = precond_event_ref.find('event').attrib['id']
            precond_event = self.tree.find("EVENT[@id='%s']" % precond_id)
            mods = self._get_mod_site(precond_event)
            if mods is None:
                warnings.warn('Skipping activity modification with missing' +\
                                'modification')
                continue

            sentence = self.get_evidence_text(event)
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id)
            self.statements.append(ActivityModification(affected_agent, mods,
                                    'increases', 'Active', evidence=ev))
            self.extracted_events['ONT::ACTIVATE'].append(event.attrib['id'])

    def get_complexes(self):
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        for event in bind_events:
            if event.attrib['id'] in self._static_events:
                continue

            sentence = self.get_evidence_text(event)
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id)

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

            self.statements.append(Complex([agent1, agent2], evidence=ev))
            self.extracted_events['ONT::BIND'].append(event.attrib['id'])

    def get_phosphorylation(self):
        phosphorylation_events = \
            self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        for event in phosphorylation_events:
            if event.attrib['id'] in self._static_events:
                continue

            enzyme = event.find(".//*[@role=':AGENT']")
            if enzyme is None:
                enzyme_agent = None
            elif enzyme.find("type").text == 'ONT::MACROMOLECULAR-COMPLEX':
                complex_id = enzyme.attrib['id']
                complex_term = self.tree.find("TERM/[@id='%s']" % complex_id)
                components = complex_term.find("components")
                terms = components.findall("term")
                term_names = []
                for t in terms:
                    term_names.append(self._get_name_by_id(t.attrib['id']))
                enzyme_name = term_names[0]
                enzyme_bound = Agent(term_names[1])
                enzyme_agent = Agent(enzyme_name,
                    bound_conditions=[BoundCondition(enzyme_bound, True)])
            else:
                enzyme_agent = self._get_agent_by_id(enzyme.attrib['id'],
                                                    event.attrib['id'])
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                warnings.warn('Skipping phosphorylation event with no '
                              'affected term.')
                continue
            affected_agent = self._get_agent_by_id(affected.attrib['id'],
                                                   event.attrib['id'])
            if affected_agent is None:
                continue
            mods = self._get_mod_site(event)

            sentence = self.get_evidence_text(event)
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id)
            # Assuming that multiple modifications can only happen in
            # distinct steps, we add a statement for each modification
            # independently

            # TODO: the first extraction here might be deprecated
            mod_types = event.findall('predicate/mods/mod/type')
            mod_types += event.findall('mods/mod/type')
            # Transphosphorylation
            if 'ONT::ACROSS' in [mt.text for mt in mod_types]:
                agent_bound = Agent(affected_agent.name)
                enzyme_agent.bound_conditions = \
                                           [BoundCondition(agent_bound, True)]
                for m in mods:
                    self.statements.append(Transphosphorylation(enzyme_agent,
                                            m.residue, m.position,
                                            evidence=ev))
            # Dephosphorylation
            elif 'ONT::MANNER-UNDO' in [mt.text for mt in mod_types]:
                for m in mods:
                    self.statements.append(Dephosphorylation(enzyme_agent,
                                            affected_agent,
                                            m.residue, m.position,
                                            evidence=ev))
            # Autophosphorylation
            elif enzyme_agent is not None and\
                (enzyme.attrib['id'] == affected.attrib['id']):
                for m in mods:
                    self.statements.append(Autophosphorylation(enzyme_agent,
                                            m.residue, m.position,
                                            evidence=ev))
            # Regular phosphorylation
            else:
                if mods is None:
                    continue
                for m in mods:
                    self.statements.append(Phosphorylation(enzyme_agent,
                                            affected_agent,
                                            m.residue, m.position,
                                            evidence=ev))
            self.extracted_events['ONT::PHOSPHORYLATION'].append(
                                                            event.attrib['id'])

    def _get_agent_by_id(self, entity_id, event_id):
        term = self.tree.find("TERM/[@id='%s']" % entity_id)
        if term is None:
            return None

        # Extract database references
        dbid = term.attrib.get('dbid')
        if dbid is None:
            db_refs_dict = {}
        else:
            dbids = dbid.split('|')
            db_refs_dict = dict([d.split(':') for d in dbids])

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
            if agent_name is None:
                return None
            agent = Agent(agent_name, db_refs=db_refs_dict)
            precond_event_ref = \
                self.tree.find("TERM/[@id='%s']/features/inevent" % entity_id)
            # Extract preconditions of the agent
            if precond_event_ref is not None:
                # Find the event describing the precondition
                preconds = precond_event_ref.findall('event')
                for precond in preconds:
                    precond_id = precond.attrib['id']
                    if precond_id == event_id:
                        warnings.warn('Circular reference to event %s.' %
                                       precond_id)
                    precond_event = self.tree.find("EVENT[@id='%s']" % 
                                                    precond_id)
                    if precond_event is None:
                        # Sometimes, if there are multiple preconditions
                        # they are numbered with <id>.1, <id>.2, etc.
                        p = self.tree.find("EVENT[@id='%s.1']" % precond_id)
                        if p is not None:
                            self.add_condition(agent, p, term)
                        p = self.tree.find("EVENT[@id='%s.2']" % precond_id)
                        if p is not None:
                            self.add_condition(agent, p, term)
                    else:
                        self.add_condition(agent, precond_event, term)
            mutations = term.findall('features/mutation')
            for mut in mutations:
                mut_id = mut.attrib.get('id')
                if mut_id is None:
                    continue
                mut_term = self.tree.find("TERM/[@id='%s']" %\
                    mut.attrib.get('id'))
                if mut_term is None:
                    continue
                mut_values = self._get_mutation(mut_term)
                if mut_values is None:
                    continue
                mc = MutCondition(mut_values[0], mut_values[1], mut_values[2])
                agent.mutations.append(mc)
        return agent

    def add_condition(self, agent, precond_event, agent_term):
        precond_event_type = precond_event.find('type').text
        # Binding precondition
        if precond_event_type == 'ONT::BIND':
            arg1 = precond_event.find('arg1')
            arg2 = precond_event.find('arg2')
            mod = precond_event.findall('mods/mod')
            bound_to_term_id = None
            if arg1 is None:
                bound_to_term_id = arg2.attrib['id']
            elif arg2 is None:
                bound_to_term_id = arg1.attrib['id']
            else:
                if arg1.attrib['id'] == agent_term.attrib['id']:
                    bound_to_term_id = arg2.attrib['id']
                else:
                    bound_to_term_id = arg1.attrib['id']

            bound_agents = []
            bound_to_term = self.tree.find("TERM/[@id='%s']" % bound_to_term_id)
            if bound_to_term.find('type').text == 'ONT::MOLECULAR-PART':
                components = bound_to_term.findall('components/component')
                for c in components:
                    bound_agent = Agent(self._get_name_by_id(c.attrib['id']))
                    bound_agents.append(bound_agent)
            else:
                bound_agents = [Agent(self._get_name_by_id(bound_to_term_id))]

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
                       agent_term.find('mods/mod[type="ONT::NEG"]')
            negation_sign = precond_event.find('predicate/negation')
            if negation_sign is not None:
                if negation_sign.text == '+':
                    neg_flag = True
            for ba in bound_agents:
                if neg_flag:
                    bc = BoundCondition(ba, False)
                else:
                    bc = BoundCondition(ba, True)
                agent.bound_conditions.append(bc)

        # Phosphorylation precondition
        elif precond_event_type == 'ONT::PHOSPHORYLATION':
            mods = self._get_mod_site(precond_event)
            agent.mods = mods

    def _find_in_term(self, term_id, path):
        tag = self.tree.find("TERM[@id='%s']/%s" % (term_id, path))
        return tag

    @staticmethod
    def _get_text(element):
        text_tag = element.find("text")
        if text_tag is None:
            return None
        text = text_tag.text
        return text

    @staticmethod
    def _get_hgnc_name(hgnc_id):
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        return hgnc_name

    @staticmethod
    def _get_valid_name(name):
        name = name.replace('-', '_')
        name = str(name.encode('utf-8').decode('ascii', 'ignore'))
        return name

    def _get_name_by_id(self, entity_id):
        entity_term = self.tree.find("TERM/[@id='%s']" % entity_id)
        if entity_term is None:
            warnings.warn('Term %s for entity not found' % entity_id)
            return None
        name = entity_term.find("name")
        if name is None:
            warnings.warn('Entity without a name')
            return None
        try:
            dbid = entity_term.attrib["dbid"]
        except:
            #warnings.warn('No grounding information for %s' % name.text)
            return self._get_valid_name(name.text)

        dbids = dbid.split('|')
        hgnc_ids = [i for i in dbids if i.startswith('HGNC')]
        up_ids = [i for i in dbids if i.startswith('UP')]


        #TODO: handle protein families like 14-3-3 with IDs like
        # XFAM:PF00244.15, FA:00007
        if hgnc_ids:
            if len(hgnc_ids) > 1:
                lisp_str = entity_term.attrib.get('lisp')
                if lisp_str is None:
                    hgnc_id = re.match(r'HGNC\:([0-9]*)',
                                       hgnc_ids[0]).groups()[0]
                else:
                    parts = lisp_str.split('(TERM :ID ')
                    scores = {}
                    for p in parts:
                        res = re.findall('HGNC::\|(.*)\|', p)
                        if res:
                            hgnc_id = res[0]
                            score = re.findall(':SCORE ([^ ]+)', p)[0]
                            scores[hgnc_id] = float(score)
                    if scores:
                        sorted_ids = sorted(scores.items(), key=operator.itemgetter(1))
                        hgnc_id = sorted_ids[-1][0]
            else:
                hgnc_id = re.match(r'HGNC\:([0-9]*)', hgnc_ids[0]).groups()[0]
            hgnc_name = self._get_hgnc_name(hgnc_id)
            return self._get_valid_name(hgnc_name)
        elif up_ids:
            if len(hgnc_ids) > 1:
                warnings.warn('%d UniProt IDs reported.' % len(up_ids))
            up_id = re.match(r'UP\:([A-Z0-9]*)', up_ids[0]).groups()[0]
            # First try to get HGNC name
            hgnc_name = up_client.get_hgnc_name(up_id)
            if hgnc_name is not None:
                return self._get_valid_name(hgnc_name)
            # Next, try to get the gene name
            gene_name = up_client.get_gene_name(up_id)
            if gene_name is not None:
                return self._get_valid_name(gene_name)
        # By default, return the text of the name tag
        name_txt = name.text.strip('|')
        return self._get_valid_name(name_txt)

    # Get all the sites recursively based on a term id.
    def _get_site_by_id(self, site_id):
        all_residues = []
        all_pos = []
        site_term = self.tree.find("TERM/[@id='%s']" % site_id)
        if site_term is None:
            # Missing site term
            return None, None

        # TODO: the 'aggregate' tag here  might be deprecated
        components = site_term.find('aggregate')
        if components is None:
            components = site_term.find('components')
        if components is not None:
            for member in components.getchildren():
                residue, pos = self._get_site_by_id(member.attrib['id'])
                all_residues += residue
                all_pos += pos
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
            elif site_type == 'ONT::MOLECULAR-DOMAIN':
                print 'Molecular domains not handled yet.'
                return None, None
            else:
                print 'Unhandled site type: %s' % site_type
                return None, None

            return (residue, ), (pos, )
        return all_residues, all_pos

    def _get_mod_site(self, event):
        # Find the modification type
        mod_type = event.find('type')
        mod_txt = mod_type.text.split('::')[1]
        mod_type_name = mod_names.get(mod_txt)
        # If it is an unknown modification type
        if mod_type_name is None:
            return None

        # Check if the event is negated
        neg = event.find('negation')
        if neg is not None and neg.text == '+':
            is_modified = False
        else:
            is_modified = True

        # Find the site of the modification
        site_tag = event.find("site")
        # If there is not site specified
        if site_tag is None:
            mc = ModCondition(mod_type_name, is_modified=is_modified)
            return [mc]
        site_id = site_tag.attrib['id']
        # Find the site TERM and get the specific residues and
        # positions
        residues, mod_pos = self._get_site_by_id(site_id)
        # If residue is missing
        if residues is None:
            mc = ModCondition(mod_type_name, is_modified=is_modified)
            return [mc]

        # Collect mods in a list
        mods = []
        for r, p in zip(residues, mod_pos):
            residue_name = get_valid_residue(r)
            if residue_name is None:
                warnings.warn('Residue name %s unknown. ' % r)
                residue_name = None
            mc = ModCondition(mod_type_name, residue_name, p, is_modified)
            mods.append(mc)
        return mods

    def _get_mutation(self, term):
        mut = term.find('mutation')
        if mut is None:
            return None
        if mut.find('type').text == 'SUBSTITUTION':
            pos = mut.find('pos').text
            aa_from = mut.find('aa-from/aa/code').text
            aa_to = mut.find('aa-to/aa/code').text
            return pos, aa_from, aa_to
        else:
            return None

    def _find_static_events(self):
        inevent_tags = self.tree.findall("TERM/features/inevent/event")
        static_events = []
        for ie in inevent_tags:
            event_id = ie.attrib['id']
            if self.tree.find("EVENT[@id='%s']" % event_id) is not None:
                static_events.append(event_id)
            else:
                # Check for events that have numbering <id>.1, <id>.2, etc.
                if self.tree.find("EVENT[@id='%s.1']" % event_id) is not None:
                    static_events.append(event_id + '.1')
                if self.tree.find("EVENT[@id='%s.2']" % event_id) is not None:
                    static_events.append(event_id + '.2')

        return static_events
