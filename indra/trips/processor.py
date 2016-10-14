from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import operator
import itertools
import collections
import xml.etree.ElementTree as ET
from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('trips')

mod_names = {
    'PHOSPHORYLATION': 'phosphorylation'
    }

molecule_types = ['ONT::GENE-PROTEIN', 'ONT::CHEMICAL', 'ONT::MOLECULE',
                 'ONT::PROTEIN', 'ONT::PROTEIN-FAMILY', 'ONT::GENE']

class TripsProcessor(object):
    """The TripsProcessor extracts INDRA Statements from a TRIPS XML.

    For more details on the TRIPS EKB XML format, see
    http://trips.ihmc.us/parser/cgi/drum

    Parameters
    ----------
    xml_string : str
        A TRIPS extraction knowledge base (EKB) in XML format as a string.

    Attributes
    ----------
    tree : xml.etree.ElementTree.Element
        An ElementTree object representation of the TRIPS EKB XML.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted from the EKB.
    doc_id : str
        The PubMed ID of the paper that the extractions are from.
    sentences : dict[str: str]
        The list of all sentences in the EKB with their IDs
    paragraphs : dict[str: str]
        The list of all paragraphs in the EKB with their IDs
    par_to_sec : dict[str: str]
        A map from paragraph IDs to their associated section types
    extracted_events : list[xml.etree.ElementTree.Element]
        A list of Event elements that have been extracted as INDRA
        Statements.
    """
    def __init__(self, xml_string):
        try:
            self.tree = ET.XML(xml_string, parser=UTB())
        except ET.ParseError:
            logger.error('Could not parse XML string')
            self.tree = None
            return
        # Get the document ID from the EKB tag. This is the PMC ID when
        # available.
        self.doc_id = self.tree.attrib.get('id')
        # Store all paragraphs and store all sentences in a data structure
        paragraph_tags = self.tree.findall('input/paragraphs/paragraph')
        sentence_tags = self.tree.findall('input/sentences/sentence')
        self.paragraphs = {p.attrib['id']: p.text for p in paragraph_tags}
        self.sentences = {s.attrib['id']: s.text for s in sentence_tags}
        self.par_to_sec = {p.attrib['id']: p.attrib.get('sec-type')
                           for p in paragraph_tags}

        self.statements = []
        self._static_events = self._find_static_events()
        self.get_all_events()
        self.extracted_events = {k:[] for k in self.all_events.keys()}
        logger.debug('All events by type')
        logger.debug('------------------')
        for k, v in self.all_events.items():
            logger.debug('%s %s' % (k, len(v)))
        logger.debug('------------------')

    def get_all_events(self):
        """Make a list of all events in the TRIPS EKB.

        The events are stored in self.all_events.
        """
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

    def get_activations(self):
        """Extract direct Activation INDRA Statements."""
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        inact_events = self.tree.findall("EVENT/[type='ONT::DEACTIVATE']")
        inact_events += self.tree.findall("EVENT/[type='ONT::INHIBIT']")
        for event in (act_events + inact_events):
            event_id = event.attrib['id']
            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec}
            ev = Evidence(source_api='trips', text=sentence,
                          pmid=self.doc_id, epistemics=epi)
            location = self._get_event_location(event)

            # Get the activating agent in the event
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                continue
            agent_id = agent.attrib.get('id')
            if agent_id is None:
                logger.debug(
                    'Skipping activation with missing activator agent')
                continue
            agent_name = self._get_name_by_id(agent_id)
            activator_agent = self._get_agent_by_id(agent_id, event_id)
            if activator_agent is None:
                continue

            # Get the activated agent in the event
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                logger.debug(
                    'Skipping activation with missing affected agent')
                continue
            affected_id = affected.attrib.get('id')
            if affected_id is None:
                logger.debug(
                    'Skipping activation with missing affected agent')
                continue

            affected_name = self._get_name_by_id(affected_id)
            if affected_name is None:
                logger.debug(
                    'Skipping activation with missing affected agent')
                continue

            affected_agent = self._get_agent_by_id(affected_id, event_id)
            if affected_agent is None:
                logger.debug(
                    'Skipping activation with missing affected agent')
                continue

            if event.find('type').text == 'ONT::ACTIVATE':
                is_activation = True
                activator_act = 'activity'
                self.extracted_events['ONT::ACTIVATE'].append(event.attrib['id'])
            elif event.find('type').text == 'ONT::INHIBIT':
                is_activation = False
                activator_act = None
                self.extracted_events['ONT::INHIBIT'].append(event.attrib['id'])
            elif event.find('type').text == 'ONT::DEACTIVATE':
                is_activation = False
                activator_act = 'activity'
                self.extracted_events['ONT::DEACTIVATE'].append(event.attrib['id'])

            for a1, a2 in _agent_list_product((activator_agent,
                                               affected_agent)):
                st = Activation(a1, activator_act, a2, 'activity',
                                is_activation=is_activation, evidence=ev)
                _stmt_location_to_agents(st, location)
                self.statements.append(st)

    def get_activations_causal(self):
        """Extract causal Activation INDRA Statements."""
        # Search for causal connectives of type ONT::CAUSE
        ccs = self.tree.findall("CC/[type='ONT::CAUSE']")
        for cc in ccs:
            factor = cc.find("arg/[@role=':FACTOR']")
            outcome = cc.find("arg/[@role=':OUTCOME']")
            # If either the factor or the outcome is missing, skip
            if factor is None or outcome is None:
                continue
            factor_id = factor.attrib.get('id')
            # Here, implicitly, we require that the factor is a TERM
            # and not an EVENT
            factor_term = self.tree.find("TERM/[@id='%s']" % factor_id)
            outcome_id = outcome.attrib.get('id')
            # Here it is implicit that the outcome is an event not
            # a TERM
            outcome_event = self.tree.find("EVENT/[@id='%s']" % outcome_id)
            if factor_term is None or outcome_event is None:
                continue
            factor_term_type = factor_term.find('type')
            # The factor term must be a molecular entity
            if factor_term_type is None or \
                factor_term_type.text not in molecule_types:
                continue
            factor_agent = self._get_agent_by_id(factor_id, None)
            if factor_agent is None:
                continue
            outcome_event_type = outcome_event.find('type')
            if outcome_event_type is None:
                continue
            # Construct evidence
            sentence = self._get_evidence_text(cc)
            sec = self._get_section(cc)
            epi = {'section_type': sec, 'direct': False}
            ev = Evidence(source_api='trips', text=sentence,
                          pmid=self.doc_id, epistemics=epi)
            location = self._get_event_location(outcome_event)
            if outcome_event_type.text == 'ONT::ACTIVATE':
                affected = outcome_event.find(".//*[@role=':AFFECTED']")
                if affected is None:
                    continue
                outcome_agent = self._get_agent_by_id(affected.attrib['id'],
                                                      outcome_id)
                if outcome_agent is None:
                    continue
                for a1, a2 in _agent_list_product((factor_agent,
                                                   outcome_agent)):
                    st = Activation(a1, 'activity',
                                    a2, 'activity', is_activation=True,
                                    evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)
            elif outcome_event_type.text == 'ONT::ACTIVITY':
                agent_tag = outcome_event.find(".//*[@role=':AGENT']")
                if agent_tag is None:
                    continue
                outcome_agent = self._get_agent_by_id(agent_tag.attrib['id'],
                                                      outcome_id)
                if outcome_agent is None:
                    continue
                for a1, a2 in _agent_list_product((factor_agent,
                                                   outcome_agent)):
                    st = Activation(a1, 'activity',
                                    a2, 'activity', is_activation=True,
                                    evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)

    def get_activations_stimulate(self):
        """Extract Activation INDRA Statements via stimulation."""
        # Search for stimulation event
        stim_events = self.tree.findall("EVENT/[type='ONT::STIMULATE']")
        for event in stim_events:
            controller = event.find("arg1/[@role=':AGENT']")
            affected = event.find("arg2/[@role=':AFFECTED']")
            # If either the controller or the affected is missing, skip
            if controller is None or affected is None:
                continue
            controller_id = controller.attrib.get('id')
            # Here, implicitly, we require that the controller is a TERM
            # and not an EVENT
            controller_term = self.tree.find("TERM/[@id='%s']" % controller_id)
            affected_id = affected.attrib.get('id')
            # Here it is implicit that the affected is an event not
            # a TERM
            affected_event = self.tree.find("EVENT/[@id='%s']" % affected_id)
            if controller_term is None or affected_event is None:
                continue
            controller_term_type = controller_term.find('type')
            # The controller term must be a molecular entity
            if controller_term_type is None or \
                controller_term_type.text not in molecule_types:
                continue
            controller_agent = self._get_agent_by_id(controller_id, None)
            if controller_agent is None:
                continue
            affected_event_type = affected_event.find('type')
            if affected_event_type is None:
                continue
            # Construct evidence
            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec, 'direct': False}
            ev = Evidence(source_api='trips', text=sentence,
                          pmid=self.doc_id, epistemics=epi)
            location = self._get_event_location(affected_event)
            if affected_event_type.text == 'ONT::ACTIVATE':
                affected = affected_event.find(".//*[@role=':AFFECTED']")
                if affected is None:
                    continue
                affected_agent = self._get_agent_by_id(affected.attrib['id'],
                                                      affected_id)
                if affected_agent is None:
                    continue
                for a1, a2 in _agent_list_product((controller_agent,
                                                   affected_agent)):
                    st = Activation(a1, 'activity',
                                    a2, 'activity', is_activation=True,
                                    evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)
            elif affected_event_type.text == 'ONT::ACTIVITY':
                agent_tag = affected_event.find(".//*[@role=':AGENT']")
                if agent_tag is None:
                    continue
                affected_agent = self._get_agent_by_id(agent_tag.attrib['id'],
                                                      affected_id)
                if affected_agent is None:
                    continue
                for a1, a2 in _agent_list_product((controller_agent,
                                                   affected_agent)):
                    st = Activation(a1, 'activity',
                                    a2, 'activity', is_activation=True,
                                    evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)

    def get_activating_mods(self):
        """Extract ActiveForm INDRA Statements based on modifications."""
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        for event in act_events:
            if event.attrib['id'] in self._static_events:
                continue
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                msg = 'Skipping activation event with no affected term.'
                logger.debug(msg)
                continue

            affected_id = affected.attrib.get('id')
            if affected_id is None:
                logger.debug(
                    'Skipping activating modification with missing' +\
                    'affected agent')
                continue

            affected_name = self._get_name_by_id(affected_id)
            if affected_name is None:
                logger.debug(
                    'Skipping activating modification with missing' +\
                    'affected agent')
                continue
            affected_agent = Agent(affected_name)
            precond_ids = self._get_precond_event_ids(affected_id)
            if not precond_ids:
                # This means that it is not an activating modification
                continue
            precond_id = precond_ids[0]
            precond_event = self.tree.find("EVENT[@id='%s']" % precond_id)
            if precond_event is None:
                continue
            mods = self._get_mod_site(precond_event)
            if mods is None:
                logger.debug('Skipping activity modification with missing' +\
                                'modification')
                continue
            affected_agent.mods = mods
            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec}
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id,
                          epistemics=epi)
            location = self._get_event_location(event)
            st = ActiveForm(affected_agent, 'activity', True, evidence=ev)
            _stmt_location_to_agents(st, location)
            self.statements.append(st)
            self.extracted_events['ONT::ACTIVATE'].append(event.attrib['id'])

    def get_complexes(self):
        """Extract Complex INDRA Statements."""
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        for event in bind_events:
            if event.attrib['id'] in self._static_events:
                continue

            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec}
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id,
                          epistemics=epi)
            location = self._get_event_location(event)

            arg1 = event.find("arg1")
            if arg1 is None or arg1.attrib.get('id') is None:
                msg = 'Skipping complex missing arg1.'
                logger.debug(msg)
                continue
            agent1 = self._get_agent_by_id(arg1.attrib['id'], event.attrib['id'])

            arg2 = event.find("arg2")
            if arg2 is None or arg2.attrib.get('id') is None:
                msg = 'Skipping complex missing arg2.'
                logger.debug(msg)
                continue
            agent2 = self._get_agent_by_id(arg2.attrib['id'], event.attrib['id'])

            # Information on binding site is either attached to the agent term
            # in a features/site tag or attached to the event itself in 
            # a site tag
            '''
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
            '''
            if agent1 is None or agent2 is None:
                logger.debug('Complex with missing members')
                continue

            for a1, a2 in _agent_list_product((agent1, agent2)):
                st = Complex([a1, a2], evidence=ev)
                _stmt_location_to_agents(st, location)
                self.statements.append(st)
            self.extracted_events['ONT::BIND'].append(event.attrib['id'])

    def get_phosphorylation(self):
        """Extract all types of phosphorylation INDRA Statements.

        The types of statements extracted here are: Phosphorylation,
        Dephosphorylation, Transphosphorylation, Autophosphorylation.
        """
        phosphorylation_events = \
            self.tree.findall("EVENT/[type='ONT::PHOSPHORYLATION']")
        for event in phosphorylation_events:
            event_id = event.attrib['id']
            if event_id in self._static_events:
                continue

            enzyme = event.find(".//*[@role=':AGENT']")
            if enzyme is None:
                enzyme_agent = None
            else:
                enzyme_id = enzyme.attrib.get('id')
                if enzyme_id is None:
                    continue
                enzyme_term = self.tree.find("TERM/[@id='%s']" % enzyme_id)
                if enzyme_term is None:
                    enzyme_agent = None
                else:
                    enzyme_agent = self._get_agent_by_id(enzyme_id, event_id)
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                logger.debug('Skipping phosphorylation event with no '
                              'affected term.')
                continue
            affected_id = affected.attrib.get('id')
            if affected_id is None:
                continue
            affected_agent = self._get_agent_by_id(affected_id, event_id)
            if affected_agent is None:
                continue
            mods = self._get_mod_site(event)

            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec}
            ev = Evidence(source_api='trips', text=sentence, pmid=self.doc_id,
                          epistemics=epi)
            location = self._get_event_location(event)
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
                    st = Transphosphorylation(enzyme_agent, m.residue,
                                              m.position, evidence=ev)
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)
            # Dephosphorylation
            elif 'ONT::MANNER-UNDO' in [mt.text for mt in mod_types]:
                for ea, aa in _agent_list_product((enzyme_agent,
                                                   affected_agent)):
                    if aa is None:
                        continue
                    for m in mods:
                        st = Dephosphorylation(ea, aa, m.residue, m.position,
                                               evidence=ev)
                        _stmt_location_to_agents(st, location)
                        self.statements.append(st)
            # Autophosphorylation
            elif enzyme_agent is not None and (enzyme_id == affected_id):
                for m in mods:
                    if isinstance(enzyme_agent, list):
                        for ea in enzyme_agent:
                            st = Autophosphorylation(ea,
                                                 m.residue, m.position,
                                                 evidence=ev)
                            _stmt_location_to_agents(st, location)
                            self.statements.append(st)
                    else:
                        st = Autophosphorylation(enzyme_agent,
                                                 m.residue, m.position,
                                                 evidence=ev)
                        _stmt_location_to_agents(st, location)
                        self.statements.append(st)
            elif affected_agent is not None and \
                'ONT::MANNER-REFL' in [mt.text for mt in mod_types]:
                for m in mods:
                    if isinstance(affected_agent, list):
                        for aa in affected_agent:
                            st = Autophosphorylation(aa,
                                                     m.residue, m.position,
                                                     evidence=ev)
                            _stmt_location_to_agents(st, location)
                            self.statements.append(st)
                    else:
                        st = Autophosphorylation(affected_agent,
                                                 m.residue, m.position,
                                                 evidence=ev)
                        _stmt_location_to_agents(st, location)
                        self.statements.append(st)

            # Regular phosphorylation
            else:
                if mods is None:
                    continue
                for ea, aa in _agent_list_product((enzyme_agent,
                                                   affected_agent)):
                    if aa is None:
                        continue
                    for m in mods:
                        st = Phosphorylation(ea, aa, m.residue, m.position,
                                             evidence=ev)
                        _stmt_location_to_agents(st, location)
                        self.statements.append(st)
            self.extracted_events['ONT::PHOSPHORYLATION'].append(
                                                            event.attrib['id'])

    def get_translocation(self):
        translocation_events = \
            self.tree.findall("EVENT/[type='ONT::TRANSLOCATE']")
        for event in translocation_events:
            event_id = event.attrib['id']
            # Get Agent which translocates
            agent_tag = event.find(".//*[@role=':AGENT']")
            if agent_tag is None:
                continue
            agent_id = agent_tag.attrib.get('id')
            agent = self._get_agent_by_id(agent_id, event_id)
            if agent is None:
                continue
            # Get from location
            from_loc_tag = event.find("from-location")
            if from_loc_tag is None:
                from_location = None
            else:
                from_loc_id = from_loc_tag.attrib.get('id')
                from_location = self._get_cell_loc_by_id(from_loc_id)
            # Get to location
            to_loc_tag = event.find("to-location")
            if to_loc_tag is None:
                to_location = None
            else:
                to_loc_id = to_loc_tag.attrib.get('id')
                to_location = self._get_cell_loc_by_id(to_loc_id)
            # Get evidence
            sentence = self._get_evidence_text(event)
            sec = self._get_section(event)
            epi = {'section_type': sec}
            ev = Evidence(source_api='trips', text=sentence,
                          pmid=self.doc_id, epistemics=epi)
            if isinstance(agent, list):
                for aa in agent:
                    st = Translocation(aa, from_location,
                                       to_location, evidence=ev)
                    self.statements.append(st)
            else:
                st = Translocation(agent, from_location,
                                   to_location, evidence=ev)
                self.statements.append(st)

    def _get_cell_loc_by_id(self, term_id):
        term = self.tree.find("TERM/[@id='%s']" % term_id)
        if term is None:
            return None
        term_type = term.find("type").text
        name = term.find("name")
        if name is None:
            return None
        else:
            name = name.text
        if term_type != 'ONT::CELL-PART':
            return None
        # If it is a cellular location, try to look up and return
        # the standard name from GO
        dbid = term.attrib.get('dbid')
        dbids = dbid.split('|')
        db_refs_dict = dict([d.split(':') for d in dbids])
        goid = db_refs_dict.get('GO')
        if goid is not None:
            try:
                loc_name = get_valid_location('GO:' + goid)
                return loc_name
            except InvalidLocationError:
                pass
        # Try to get the same from UP
        upid = db_refs_dict.get('UP')
        if upid is not None and upid.startswith('SL'):
            loc_name = up_client.uniprot_subcell_loc.get(upid)
            if loc_name is not None:
                try:
                    loc_name = get_valid_location(loc_name.lower())
                    return loc_name
                except InvalidLocationError:
                    pass
        # Check if the raw name is a valid cellular component
        if name is not None:
            try:
                loc_name = get_valid_location(name.lower())
                return loc_name
            except InvalidLocationError:
                pass
        msg = 'Location %s is not a valid GO cellular component' % name
        logger.debug(msg)
        return None

    def _get_event_location(self, event_term):
        location = event_term.find('location')
        if location is None:
            return None
        loc_id = location.get('id')
        loc = self._get_cell_loc_by_id(loc_id)
        return loc

    def _get_agent_by_id(self, entity_id, event_id):
        term = self.tree.find("TERM/[@id='%s']" % entity_id)
        if term is None:
            return None

        # Check if the term is an aggregate
        members = term.findall('aggregate/member')
        if members:
            op = term.find('aggregate').attrib.get('operator')
            if op != 'AND':
                logger.debug('Skipping aggregate with operator %s' % op)
                return None
            member_ids = [m.attrib.get('id') for m in members]
            member_agents = [self._get_agent_by_id(m, event_id)
                             for m in member_ids]
            return member_agents

        db_refs_dict = self._get_db_refs(term)

        agent_text_tag = term.find('name')
        if agent_text_tag is not None:
            agent_text = agent_text_tag.text
            db_refs_dict['TEXT'] = agent_text

        # If the entity is a complex
        if term.find("type").text == 'ONT::MACROMOLECULAR-COMPLEX':
            complex_id = entity_id
            complex_term = self.tree.find("TERM/[@id='%s']" % complex_id)
            components = complex_term.find("components")
            if components is None:
                logger.debug('Complex without components')
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
            precond_ids = self._get_precond_event_ids(entity_id)
            if precond_ids:
                for precond_id in precond_ids:
                    if precond_id == event_id:
                        logger.debug('Circular reference to event %s.' %
                                       precond_id)
                    precond_event = self.tree.find("EVENT[@id='%s']" % 
                                                    precond_id)
                    if precond_event is None:
                        # Sometimes, if there are multiple preconditions
                        # they are numbered with <id>.1, <id>.2, etc.
                        p = self.tree.find("EVENT[@id='%s.1']" % precond_id)
                        if p is not None:
                            self._add_condition(agent, p, term)
                        p = self.tree.find("EVENT[@id='%s.2']" % precond_id)
                        if p is not None:
                            self._add_condition(agent, p, term)
                    else:
                        self._add_condition(agent, precond_event, term)
            # Get mutations
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
                try:
                    mc = MutCondition(mut_values[0], mut_values[1], mut_values[2])
                except InvalidResidueError:
                    logger.error('Invalid residue in mutation condition.')
                    continue
                agent.mutations.append(mc)
        # Get location
        location = term.find('features/location')
        if location is not None:
            loc_id = location.attrib.get('id')
            loc = self._get_cell_loc_by_id(loc_id)
            agent.location = loc
        # Get activity
        activity = term.find('features/active')
        if activity is not None:
            if activity.text.lower() == 'true':
                agent.active = 'activity'

        return agent

    @staticmethod
    def _get_db_refs(term):
        # Extract database references
        dbid = term.attrib.get('dbid')
        if dbid is None:
            db_refs_dict = {}
            if term.find('type').text == 'ONT::PROTEIN-FAMILY':
                members = term.findall('members/member')
                dbids = []
                for m in members:
                    dbid = m.attrib.get('dbid')
                    parts = dbid.split(':')
                    dbids.append({parts[0]: parts[1]})
                db_refs_dict = {'PFAM-DEF': dbids}
        else:
            drum_terms = term.findall('drum-terms/drum-term')
            if drum_terms:
                scores = {}
                score_started = False
                for dt in drum_terms:
                    dbid_str = dt.attrib.get('dbid')
                    match_score = dt.attrib.get('match-score')
                    if not score_started:
                        if match_score is not None:
                            score_started = True
                        else:
                            # This is a match before other scored terms so we
                            # default to 1.0
                            match_score = 1.0
                    else:
                        if match_score is None:
                            # This is a match after other scored matches
                            # default to a small value
                            match_score = 0.1

                    if dbid_str is None:
                        db_refs_dict = {}
                        if term.find('type').text == 'ONT::PROTEIN-FAMILY':
                            members = term.findall('members/member')
                            dbids = []
                            for m in members:
                                dbid = m.attrib.get('dbid')
                                dbids.append(dbid)
                            key_name = 'PFAM-DEF:' + '|'.join(dbids)
                            scores[key_name] = float(match_score)
                    else:
                        scores[dbid_str] = float(match_score)
                    xr_tags = dt.findall('xrefs/xref')
                    for xrt in xr_tags:
                        dbid_str = xrt.attrib.get('dbid')
                        scores[dbid_str] = float(match_score)
                sorted_db_refs = sorted(scores.items(),
                                        key=operator.itemgetter(1),
                                        reverse=True)
                db_refs_dict = {}
                for dbid_str, _ in sorted_db_refs:
                    dbname, dbid = dbid_str.split(':')
                    if not db_refs_dict.get(dbname):
                        if dbname == 'PFAM-DEF':
                            dbids = [{p[0]: p[1]} for p in dbid.split('|')]
                            db_refs_dict[dbname] = dbids
                        else:
                            db_refs_dict[dbname] = dbid

            else:
                dbids = dbid.split('|')
                db_refs_dict = {}
                for dbname, dbid in [d.split(':') for d in dbids]:
                    if not db_refs_dict.get(dbname):
                        db_refs_dict[dbname] = dbid
        return db_refs_dict

    def _add_condition(self, agent, precond_event, agent_term):
        precond_event_type = precond_event.find('type').text
        # Binding precondition
        if precond_event_type == 'ONT::BIND':
            arg1 = precond_event.find('arg1')
            arg2 = precond_event.find('arg2')
            mod = precond_event.findall('mods/mod')
            bound_to_term_id = None
            if arg1 is None:
                bound_to_term_id = arg2.attrib.get('id')
            elif arg2 is None:
                bound_to_term_id = arg1.attrib.get('id')
            else:
                arg1_id = arg1.attrib.get('id')
                arg2_id = arg2.attrib.get('id')
                if arg1_id == agent_term.attrib['id']:
                    bound_to_term_id = arg2_id
                else:
                    bound_to_term_id = arg1_id

            bound_agents = []
            if bound_to_term_id is not None:
                bound_to_term = self.tree.find("TERM/[@id='%s']" % bound_to_term_id)
                if bound_to_term.find('type').text == 'ONT::MOLECULAR-PART':
                    components = bound_to_term.findall('components/component')
                    for c in components:
                        bound_agent_name = self._get_name_by_id(c.attrib['id'])
                        if bound_agent_name is not None:
                            bound_agent = Agent(bound_agent_name)
                            bound_agents.append(bound_agent)
                else:
                    bound_agent_name = self._get_name_by_id(bound_to_term_id)
                    if bound_agent_name is not None:
                        bound_agents = [Agent(bound_agent_name)]

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
        name = name.replace('/', '_')
        name = str(name.encode('utf-8').decode('ascii', 'ignore'))
        return name

    def _get_name_by_id(self, entity_id):
        entity_term = self.tree.find("TERM/[@id='%s']" % entity_id)
        if entity_term is None:
            logger.debug('Term %s for entity not found' % entity_id)
            return None
        name = entity_term.find("name")
        if name is None:
            logger.debug('Entity without a name')
            return None
        db_refs = self._get_db_refs(entity_term)
        if not db_refs:
            return self._get_valid_name(name.text)

        #TODO: handle protein families like 14-3-3 with IDs like
        # XFAM:PF00244.15, FA:00007
        hgnc_id = db_refs.get('HGNC')
        up_id = db_refs.get('UP')
        if hgnc_id:
            hgnc_name = self._get_hgnc_name(hgnc_id)
            return self._get_valid_name(hgnc_name)
        elif up_id:
            # First to get gene name
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
            site_name_tag = site_term.find("name")
            if site_name_tag is not None:
                site_name = site_name_tag.text
            if site_type == 'ONT::MOLECULAR-SITE':
                residue = site_term.find('features/site/code')
                if residue is not None:
                    residue = residue.text.upper()
                pos = site_term.find('features/site/pos')
                if pos is not None:
                    pos = pos.text.upper()
            elif site_type == 'ONT::RESIDUE':
                # Example name: TYROSINE-RESIDUE
                if site_name is not None:
                    residue = site_name.split('-')[0]
                else:
                    residue = None
                pos = None
            elif site_type == 'ONT::AMINO-ACID':
                residue = site_name
                pos = None
            elif site_type == 'ONT::MOLECULAR-DOMAIN':
                logger.debug('Molecular domains not handled yet.')
                return None, None
            else:
                logger.debug('Unhandled site type: %s' % site_type)
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
            try:
                residue_name = get_valid_residue(r)
            except InvalidResidueError:
                logger.debug('Invalid residue name %s' % r)
                residue_name = None
            mc = ModCondition(mod_type_name, residue_name, p, is_modified)
            mods.append(mc)
        return mods

    def _get_mutation(self, term):
        mut = term.find('mutation')
        if mut is None or mut.find('type') is None:
            return None
        if mut.find('type').text == 'SUBSTITUTION':
            pos_tag = mut.find('pos')
            if pos_tag is not None:
                pos = pos_tag.text
            else:
                pos = None
            aa_from_tag = mut.find('aa-from/aa/code')
            if aa_from_tag is not None:
                aa_from = aa_from_tag.text
            else:
                aa_from = None
            aa_to_tag = mut.find('aa-to/aa/code')
            if aa_to_tag is not None:
                aa_to = aa_to_tag.text
            else:
                aa_to = None
            return pos, aa_from, aa_to
        else:
            return None

    def _get_evidence_text(self, event_tag):
        """Extract the evidence for an event.

        Pieces of text linked to an EVENT are fragments of a sentence. The
        EVENT refers to the paragraph ID and the "uttnum", which corresponds 
        to a sentence ID. Here we find and return the full sentence from which 
        the event was taken.
        """
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

    def _get_section(self, event_tag):
        par_id = event_tag.attrib.get('paragraph')
        sec = self.par_to_sec.get(par_id)
        return sec

    def _get_precond_event_ids(self, term_id):
        precond_ids = []
        precond_event_ref = \
            self.tree.find("TERM/[@id='%s']/features/inevent" % term_id)
        if precond_event_ref is not None:
            preconds = precond_event_ref.findall('event')
            precond_ids += [p.attrib.get('id') for p in preconds]
        precond_event_refs = \
            self.tree.findall("TERM/[@id='%s']/features/ptm" % term_id)
        precond_ids += [p.attrib.get('event') for p in precond_event_refs]
        return precond_ids

    def _find_static_events(self):
        # Find sub-EVENTs that TERMs refer to
        inevent_tags = self.tree.findall("TERM/features/inevent/event")
        ptm_tags = self.tree.findall("TERM/features/ptm")
        notptm_tags = self.tree.findall("TERM/features/not-ptm")
        sub_event_ids = [t.attrib.get('id') for t in inevent_tags]
        sub_event_ids += [t.attrib.get('event') for t in ptm_tags]
        sub_event_ids += [t.attrib.get('event') for t in notptm_tags]
        static_events = []
        for event_id in sub_event_ids:
            if event_id == 'V2260949':
                import ipdb; ipdb.set_trace()
            event_tag = self.tree.find("EVENT[@id='%s']" % event_id)
            if event_tag is not None:
                # If an affected TERM in the primary event has the same event
                # specified as a not-ptm, that doesn't count as a static
                # event. Therefore we let these events go through.
                affected = event_tag.find(".//*[@role=':AFFECTED']")
                if affected is not None:
                    affected_id = affected.attrib.get('id')
                    enp = self.tree.find("TERM[@id='%s']/not-features/ptm" %
                                         affected_id)
                    if (enp is not None and
                        enp.attrib.get('event') == event_id):
                        continue
                static_events.append(event_id)
            else:
                # Check for events that have numbering <id>.1, <id>.2, etc.
                if self.tree.find("EVENT[@id='%s.1']" % event_id) is not None:
                    static_events.append(event_id + '.1')
                if self.tree.find("EVENT[@id='%s.2']" % event_id) is not None:
                    static_events.append(event_id + '.2')

        return static_events

def _stmt_location_to_agents(stmt, location):
    """Apply an event location to the Agents in the corresponding Statement.

    If a Statement is in a given location we represent that by requiring all
    Agents in the Statement to be in that location.
    """
    agents = stmt.agent_list()
    for a in agents:
        if a is not None:
            a.location = location

def _agent_list_product(lists):
    def _listify(lst):
        if not isinstance(lst, collections.Iterable):
            return [lst]
        else:
            return lst
    ll = [_listify(l) for l in lists]
    return itertools.product(*ll)
