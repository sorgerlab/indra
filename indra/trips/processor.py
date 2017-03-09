from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import re
import logging
import operator
import itertools
import collections
import xml.etree.ElementTree as ET
from indra.util import read_unicode_csv
from indra.statements import *
import indra.databases.hgnc_client as hgnc_client
import indra.databases.uniprot_client as up_client
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('trips')

mod_names = {
    'ONT::PHOSPHORYLATION': 'phosphorylation',
    'ONT::UBIQUITINATION': 'ubiquitination',
    'ONT::RIBOSYLATION': 'ribosylation',
    'ONT::ACETYLATION': 'acetylation',
    'ONT::HYDROXYLATION': 'hydroxylation',
    'ONT::FARNESYLATION': 'farnesylation'
    }

protein_types = ['ONT::GENE-PROTEIN', 'ONT::CHEMICAL', 'ONT::MOLECULE',
                 'ONT::PROTEIN', 'ONT::PROTEIN-FAMILY', 'ONT::GENE',
                 'ONT::MACROMOLECULAR-COMPLEX']

molecule_types = protein_types + \
    ['ONT::CHEMICAL', 'ONT::MOLECULE', 'ONT::SUBSTANCE',
     'ONT::PHARMACOLOGIC-SUBSTANCE']


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
        self._isolated_terms = self._find_isolated_terms()
        self._subsumed_events = []
        self.all_events = {}
        self.get_all_events()
        self.extracted_events = {k: [] for k in self.all_events.keys()}
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
            if event_id in self._static_events:
                continue
            # Get the activating agent in the event
            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                continue
            agent_id = agent.attrib.get('id')
            if agent_id is None:
                logger.debug(
                    'Skipping activation with missing activator agent')
                continue
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

            affected_agent = self._get_agent_by_id(affected_id, event_id)
            if affected_agent is None:
                logger.debug(
                    'Skipping activation with missing affected agent')
                continue

            is_activation = True
            if _is_type(event, 'ONT::ACTIVATE'):
                self._add_extracted('ONT::ACTIVATE', event.attrib['id'])
            elif _is_type(event, 'ONT::INHIBIT'):
                is_activation = False
                self._add_extracted('ONT::INHIBIT', event.attrib['id'])
            elif _is_type(event, 'ONT::DEACTIVATE'):
                is_activation = False
                self._add_extracted('ONT::DEACTIVATE', event.attrib['id'])

            ev = self._get_evidence(event)
            location = self._get_event_location(event)

            for a1, a2 in _agent_list_product((activator_agent,
                                               affected_agent)):
                if is_activation:
                    st = Activation(a1, a2, evidence=[ev])
                else:
                    st = Inhibition(a1, a2, evidence=[ev])
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
            ev = self._get_evidence(cc)
            ev.epistemics['direct'] = False
            location = self._get_event_location(outcome_event)
            if outcome_event_type.text in ['ONT::ACTIVATE', 'ONT::ACTIVITY',
                                           'ONT::DEACTIVATE']:
                if outcome_event_type.text in ['ONT::ACTIVATE',
                                               'ONT::DEACTIVATE']:
                    agent_tag = outcome_event.find(".//*[@role=':AFFECTED']")
                elif outcome_event_type.text == 'ONT::ACTIVITY':
                    agent_tag = outcome_event.find(".//*[@role=':AGENT']")
                if agent_tag is None or agent_tag.attrib.get('id') is None:
                    continue
                outcome_agent = self._get_agent_by_id(agent_tag.attrib['id'],
                                                      outcome_id)
                if outcome_agent is None:
                    continue
                if outcome_event_type.text == 'ONT::DEACTIVATE':
                    is_activation = False
                else:
                    is_activation = True
                for a1, a2 in _agent_list_product((factor_agent,
                                                   outcome_agent)):
                    if is_activation:
                        st = Activation(a1, a2, evidence=[ev])
                    else:
                        st = Inhibition(a1, a2, evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)

    def get_activations_stimulate(self):
        """Extract Activation INDRA Statements via stimulation."""
        # TODO: extract to other patterns:
        # - Stimulation by EGF activates ERK
        # - Stimulation by EGF leads to ERK activation
        # Search for stimulation event
        stim_events = self.tree.findall("EVENT/[type='ONT::STIMULATE']")
        for event in stim_events:
            event_id = event.attrib.get('id')
            if event_id in self._static_events:
                continue
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
            ev = self._get_evidence(event)
            ev.epistemics['direct'] = False
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
                    st = Activation(a1, a2, evidence=[ev])
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
                    st = Activation(a1, a2, evidence=[ev])
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)

    def get_degradations(self):
        """Extract Degradation INDRA Statements."""
        deg_events = self.tree.findall("EVENT/[type='ONT::CONSUME']")
        for event in deg_events:
            if event.attrib['id'] in self._static_events:
                continue
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                msg = 'Skipping degradation event with no affected term.'
                logger.debug(msg)
                continue

            # Make sure the degradation is affecting a molecule type
            # Temporarily removed for CwC compatibility with no type tag
            #affected_type = affected.find('type')
            #if affected_type is None or \
            #   affected_type.text not in molecule_types:
            #    continue

            affected_id = affected.attrib.get('id')
            if affected_id is None:
                logger.debug(
                    'Skipping degradation event with missing affected agent')
                continue

            affected_agent = self._get_agent_by_id(affected_id,
                                                   event.attrib['id'])
            if affected_agent is None:
                logger.debug(
                    'Skipping degradation event with missing affected agent')
                continue

            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                agent_agent = None
            else:
                agent_id = agent.attrib.get('id')
                if agent_id is None:
                    agent_agent = None
                else:
                    agent_agent = self._get_agent_by_id(agent_id,
                                                        event.attrib['id'])

            ev = self._get_evidence(event)
            location = self._get_event_location(event)
            for subj, obj in \
                    _agent_list_product((agent_agent, affected_agent)):
                st = DecreaseAmount(subj, obj, evidence=ev)
                _stmt_location_to_agents(st, location)
                self.statements.append(st)
            self._add_extracted(_get_type(event), event.attrib['id'])

    def get_syntheses(self):
        """Extract IncreaseAmount INDRA Statements."""
        syn_events = self.tree.findall("EVENT/[type='ONT::PRODUCE']")
        syn_events += self.tree.findall("EVENT/[type='ONT::TRANSCRIBE']")
        for event in syn_events:
            if event.attrib['id'] in self._static_events:
                continue
            if event.attrib['id'] in self._subsumed_events:
                continue
            affected = event.find(".//*[@role=':AFFECTED-RESULT']")
            if affected is None:
                msg = 'Skipping synthesis event with no affected term.'
                logger.debug(msg)
                continue

            # Make sure the synthesis is affecting a molecule type
            # Temporarily removed for CwC compatibility with no type tag
            # affected_type = affected.find('type')
            # if affected_type is None or \
            #   affected_type.text not in molecule_types:
            #    continue

            affected_id = affected.attrib.get('id')
            if affected_id is None:
                logger.debug(
                    'Skipping synthesis event with missing affected agent')
                continue

            affected_agent = self._get_agent_by_id(affected_id,
                                                   event.attrib['id'])
            if affected_agent is None:
                logger.debug(
                    'Skipping synthesis event with missing affected agent')
                continue

            agent = event.find(".//*[@role=':AGENT']")
            if agent is None:
                agent_agent = None
            else:
                agent_id = agent.attrib.get('id')
                if agent_id is None:
                    agent_agent = None
                else:
                    agent_agent = self._get_agent_by_id(agent_id,
                                                        event.attrib['id'])
                    '''
                    if _get_type(event) == 'ONT::TRANSCRIBE':
                        if agent_agent is not None:
                            agent_agent.activity = \
                                    ActivityCondition('transcription', True)
                    '''
            ev = self._get_evidence(event)
            location = self._get_event_location(event)
            for subj, obj in \
                    _agent_list_product((agent_agent, affected_agent)):
                st = IncreaseAmount(subj, obj, evidence=ev)
                _stmt_location_to_agents(st, location)
                self.statements.append(st)
            self._add_extracted(_get_type(event), event.attrib['id'])

    def get_regulate_amounts(self):
        """Extract Increase/DecreaseAmount Statements."""
        pos_events = []
        neg_events = []
        pattern = "EVENT/[type='ONT::STIMULATE']/arg2/[type='ONT::TRANSCRIBE']/.."
        pos_events += self.tree.findall(pattern)
        pattern = "EVENT/[type='ONT::INCREASE']/arg2/[type='ONT::TRANSCRIBE']/.."
        pos_events += self.tree.findall(pattern)
        pattern = "EVENT/[type='ONT::INHIBIT']/arg2/[type='ONT::TRANSCRIBE']/.."
        neg_events += self.tree.findall(pattern)
        pattern = "EVENT/[type='ONT::DECREASE']/arg2/[type='ONT::TRANSCRIBE']/.."
        neg_events += self.tree.findall(pattern)
        # Look at polarity
        pattern = "EVENT/[type='ONT::MODULATE']/arg2/[type='ONT::TRANSCRIBE']/.."
        mod_events = self.tree.findall(pattern)
        for event in mod_events:
            pol = event.find('polarity')
            if pol is not None:
                if pol.text == 'ONT::POSITIVE':
                    pos_events.append(event)
                elif pol.text == 'ONT::NEGATIVE':
                    neg_events.append(event)
        combs = zip([pos_events, neg_events], [IncreaseAmount, DecreaseAmount])
        for events, cls in combs:
            for event in events:
                # The agent has to exist and be a protein type
                agent = event.find(".//*[@role=':AGENT']")
                if agent is None:
                    continue
                if agent.find('type') is None or \
                    (agent.find('type').text not in protein_types):
                    continue
                agent_id = agent.attrib.get('id')
                if agent_id is None:
                    continue
                agent_agent = self._get_agent_by_id(agent_id,
                                                    event.attrib['id'])

                # The affected, we already know is ONT::TRANSCRIPTION
                affected_arg = event.find(".//*[@role=':AFFECTED']")
                affected_id = affected_arg.attrib.get('id')
                affected_event = self.tree.find("EVENT/[@id='%s']" %
                                                affected_id)
                if affected_event is None:
                    continue
                affected = \
                    affected_event.find(".//*[@role=':AFFECTED-RESULT']")
                if affected is None:
                    affected = \
                        affected_event.find(".//*[@role=':AFFECTED']")
                    if affected is None:
                        continue
                affected_id = affected.attrib.get('id')
                if affected_id is None:
                    continue
                affected_agent = \
                        self._get_agent_by_id(affected_id,
                                              affected_event.attrib['id'])
                ev = self._get_evidence(event)
                location = self._get_event_location(event)
                for subj, obj in \
                        _agent_list_product((agent_agent, affected_agent)):
                    st = cls(subj, obj, evidence=ev)
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)
                self._add_extracted(_get_type(event), event.attrib['id'])
                self._subsumed_events.append(affected_event.attrib['id'])


    def get_active_forms(self):
        """Extract ActiveForm INDRA Statements."""
        act_events = self.tree.findall("EVENT/[type='ONT::ACTIVATE']")
        def _agent_is_basic(agent):
            if not agent.mods and not agent.mutations \
                    and not agent.bound_conditions and not agent.location:
                return True
            return False
        for event in act_events:
            if event.attrib['id'] in self._static_events:
                continue
            agent = event.find(".//*[@role=':AGENT']")
            if agent is not None:
                # In this case this is not an ActiveForm statement
                continue
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                msg = 'Skipping active form event with no affected term.'
                logger.debug(msg)
                continue

            affected_id = affected.attrib.get('id')
            if affected_id is None:
                logger.debug(
                    'Skipping active form event with missing affected agent')
                continue

            affected_agent = self._get_agent_by_id(affected_id,
                                                   event.attrib['id'])
            # If it is a list of agents, skip them for now
            if not isinstance(affected_agent, Agent):
                continue
            if _agent_is_basic(affected_agent):
                continue
            # The affected agent has to be protein-like type
            affected_type = affected.find('type')
            if affected_type is None or \
               affected_type.text not in protein_types:
                continue
            # If the Agent state is at the base state then this is not an
            # ActiveForm statement
            if _is_base_agent_state(affected_agent):
                continue
            ev = self._get_evidence(event)
            location = self._get_event_location(event)
            st = ActiveForm(affected_agent, 'activity', True, evidence=ev)
            _stmt_location_to_agents(st, location)
            self.statements.append(st)
            self._add_extracted('ONT::ACTIVATE', event.attrib['id'])

    def get_active_forms_state(self):
        """Extract ActiveForm INDRA Statements."""
        for term in self._isolated_terms:
            act = term.find('features/active')
            if act is None:
                continue
            if act.text == 'TRUE':
                is_active = True
            elif act.text == 'FALSE':
                is_active = False
            else:
                logger.warning('Unhandled term activity feature %s' % act.text)
            agent = self._get_agent_by_id(term.attrib['id'], None)
            # Skip aggregates for now
            if not isinstance(agent, Agent):
                continue
            # If the Agent state is at the base state then this is not an
            # ActiveForm statement
            if _is_base_agent_state(agent):
                continue
            # Remove the activity flag since it's irrelevant here
            agent.activity = None
            text_term = term.find('text')
            if text_term is not None:
                ev_text = text_term.text
            else:
                ev_text = None
            ev = Evidence(source_api='trips', text=ev_text, pmid=self.doc_id)
            st = ActiveForm(agent, 'activity', is_active, evidence=[ev])
            self.statements.append(st)

    def get_complexes(self):
        """Extract Complex INDRA Statements."""
        bind_events = self.tree.findall("EVENT/[type='ONT::BIND']")
        bind_events += self.tree.findall("EVENT/[type='ONT::INTERACT']")
        for event in bind_events:
            if event.attrib['id'] in self._static_events:
                continue

            arg1 = event.find("arg1")
            arg2 = event.find("arg2")
            if (arg1 is None or arg1.attrib.get('id') is None) or \
               (arg2 is None or arg2.attrib.get('id') is None):
                logger.debug('Skipping complex with less than 2 members')
                continue

            agent1 = self._get_agent_by_id(arg1.attrib['id'],
                                           event.attrib['id'])
            agent2 = self._get_agent_by_id(arg2.attrib['id'],
                                           event.attrib['id'])
            if agent1 is None or agent2 is None:
                logger.debug('Skipping complex with less than 2 members')
                continue

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
            ev = self._get_evidence(event)
            location = self._get_event_location(event)

            for a1, a2 in _agent_list_product((agent1, agent2)):
                st = Complex([a1, a2], evidence=ev)
                _stmt_location_to_agents(st, location)
                self.statements.append(st)
            self._add_extracted(_get_type(event), event.attrib['id'])

    def get_modifications(self):
        """Extract all types of Modification INDRA Statements."""
        mod_event_types = mod_names.keys()
        mod_events = []
        for mod_event_type in mod_event_types:
            events = self.tree.findall("EVENT/[type='%s']" % mod_event_type)
            mod_events += events
        for event in mod_events:
            event_id = event.attrib['id']
            event_type = _get_type(event)
            if event_id in self._static_events:
                continue

            # Get enzyme Agent
            enzyme = event.find(".//*[@role=':AGENT']")
            if enzyme is None:
                enzyme_agent = None
            else:
                enzyme_id = enzyme.attrib.get('id')
                if enzyme_id is None:
                    continue
                enzyme_agent = self._get_agent_by_id(enzyme_id, event_id)

            # Get substrate Agent
            affected = event.find(".//*[@role=':AFFECTED']")
            if affected is None:
                logger.debug('Skipping modification event with no '
                             'affected term.')
                continue
            affected_id = affected.attrib.get('id')
            if affected_id is None:
                continue
            affected_agent = self._get_agent_by_id(affected_id, event_id)
            if affected_agent is None:
                logger.debug('Skipping modification event with no '
                             'affected term.')
                continue

            # Get modification sites
            mods = self._get_modification(event)

            # Get evidence and location
            ev = self._get_evidence(event)
            location = self._get_event_location(event)

            mod_types = event.findall('mods/mod/type')

            # Trans and Auto are unique to Phosphorylation
            if _is_type(event, 'ONT::PHOSPHORYLATION'):
                # Transphosphorylation
                if enzyme_agent is not None and \
                        'ONT::ACROSS' in [mt.text for mt in mod_types]:
                    agent_bound = Agent(affected_agent.name)
                    enzyme_agent.bound_conditions = \
                        [BoundCondition(agent_bound, True)]
                    for m in mods:
                        st = Transphosphorylation(enzyme_agent, m.residue,
                                                  m.position, evidence=[ev])
                        _stmt_location_to_agents(st, location)
                        self.statements.append(st)
                    continue
                # Autophosphorylation
                elif enzyme_agent is not None and (enzyme_id == affected_id):
                    for m in mods:
                        if isinstance(enzyme_agent, list):
                            for ea in enzyme_agent:
                                st = Autophosphorylation(ea,
                                                     m.residue, m.position,
                                                     evidence=[ev])
                                _stmt_location_to_agents(st, location)
                                self.statements.append(st)
                        else:
                            st = Autophosphorylation(enzyme_agent,
                                                     m.residue, m.position,
                                                     evidence=[ev])
                            _stmt_location_to_agents(st, location)
                            self.statements.append(st)
                    continue
                elif affected_agent is not None and \
                     'ONT::MANNER-REFL' in [mt.text for mt in mod_types]:
                    for m in mods:
                        if isinstance(affected_agent, list):
                            for aa in affected_agent:
                                st = Autophosphorylation(aa,
                                                         m.residue, m.position,
                                                         evidence=[ev])
                                _stmt_location_to_agents(st, location)
                                self.statements.append(st)
                        else:
                            st = Autophosphorylation(affected_agent,
                                                     m.residue, m.position,
                                                     evidence=[ev])
                            _stmt_location_to_agents(st, location)
                            self.statements.append(st)
                    continue

            mod = mod_names.get(event_type)
            if 'ONT::MANNER-UNDO' in [mt.text for mt in mod_types]:
                mod_stmt = modclass_to_inverse[modtype_to_modclass[mod]]
            else:
                mod_stmt = modtype_to_modclass[mod]
            for ea, aa in _agent_list_product((enzyme_agent, affected_agent)):
                if aa is None:
                    continue
                for m in mods:
                    st = mod_stmt(ea, aa, m.residue, m.position, evidence=ev)
                    _stmt_location_to_agents(st, location)
                    self.statements.append(st)
            self._add_extracted(event_type, event.attrib['id'])

    def get_translocation(self):
        translocation_events = \
            self.tree.findall("EVENT/[type='ONT::TRANSLOCATE']")
        for event in translocation_events:
            event_id = event.attrib['id']
            if event_id in self._static_events:
                continue
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
            if from_location is None and to_location is None:
                continue
            # Get evidence
            ev = self._get_evidence(event)
            if isinstance(agent, list):
                for aa in agent:
                    st = Translocation(aa, from_location,
                                       to_location, evidence=ev)
                    self.statements.append(st)
            else:
                st = Translocation(agent, from_location,
                                   to_location, evidence=ev)
                self.statements.append(st)
            self._add_extracted('ONT::TRANSLOCATE', event.attrib['id'])

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
                logger.debug('Skipping aggregate with operator %s.' % op)
                return None
            member_ids = [m.attrib.get('id') for m in members]
            member_agents = []
            for member_id in member_ids:
                agent = self._get_agent_by_id(member_id, event_id)
                if agent is None:
                    logger.warning('Could not extract term %s.' %
                                   member_id)
                    continue
                if isinstance(agent, Agent):
                    member_agents.append(agent)
                else:
                    member_agents += agent
            # Handle case where the individual member extraction fails
            # to make sure we don't end up with None Agent arguments
            # in Statements
            if not member_agents:
                return None
            return member_agents

        db_refs, _ = _get_db_refs(term)

        # If the entity is a complex
        if _is_type(term, 'ONT::MACROMOLECULAR-COMPLEX'):
            components = term.findall("components/component")
            agents = []
            for component in components:
                component_id = component.attrib['id']
                agent = self._get_agent_by_id(component_id, None)
                if agent is not None:
                    agents.append(agent)
            if not agents:
                return None
            # We assume that the first agent mentioned in the description of
            # the complex is the one that mediates binding
            agent = agents[0]
            if len(agents) > 1:
                agent.bound_conditions = \
                    [BoundCondition(ag, True) for ag in agents[1:]]
        # If the entity is not a complex
        else:
            # Determine the agent name
            hgnc_id = db_refs.get('HGNC')
            up_id = db_refs.get('UP')
            be_id = db_refs.get('BE')
            agent_name = None
            # HGNC name takes precedence
            if hgnc_id:
                hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
                agent_name = hgnc_name
            # If no HGNC name (for instance non-human protein) then
            # look at UP and try to get gene name
            elif up_id:
                gene_name = up_client.get_gene_name(up_id)
                if gene_name:
                    agent_name = gene_name
            # If it is mapped to Bioentities then we standardize its name
            # to the Bioentities entry name
            elif be_id:
                agent_name = be_id
            # Otherwise, take the name of the term as agent name
            else:
                name = term.find("name")
                if name is not None:
                    agent_name = name.text
            # If after all of this, the agent name is still None
            # then we don't extract this term as an agent
            if agent_name is None:
                return None
            agent = Agent(agent_name, db_refs=db_refs)

        # Look for precondition events and apply them to the Agent
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
            mut_term = self.tree.find("TERM/[@id='%s']" % mut.attrib.get('id'))
            if mut_term is None:
                continue
            mut_values = self._get_mutation(mut_term)
            if mut_values is None:
                continue
            try:
                mc = MutCondition(mut_values[0], mut_values[1],
                                  mut_values[2])
            except InvalidResidueError:
                residues_str = '%s/%s' % (mut_values[1], mut_values[2])
                logger.error('Invalid residue in mutation condition: %s' % \
                             residues_str)
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
                agent.activity = ActivityCondition('activity', True)
            if activity.text.lower() == 'false':
                agent.activity = ActivityCondition('activity', False)
        return agent

    def _add_condition(self, agent, precond_event, agent_term):
        precond_event_type = _get_type(precond_event)

        # Modification precondition
        if precond_event_type in mod_names.keys():
            mods = self._get_modification(precond_event)
            agent.mods = mods
            return

        # Binding precondition
        if precond_event_type == 'ONT::BIND':
            arg1 = precond_event.find('arg1')
            arg2 = precond_event.find('arg2')
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
            if bound_to_term_id == agent_term.attrib['id']:
                return

            bound_agents = []
            if bound_to_term_id is not None:
                bound_to_term = self.tree.find("TERM/[@id='%s']" %
                                               bound_to_term_id)
                if bound_to_term is None:
                    pass
                elif _is_type(bound_to_term, 'ONT::CELL-PART'):
                    # We currently don't handle binding to cellular components
                    # TODO: possibly handle this as location
                    pass
                elif _is_type(bound_to_term, 'ONT::MOLECULAR-PART'):
                    components = bound_to_term.findall('components/component')
                    for c in components:
                        bound_agent = \
                            self._get_basic_agent_by_id(c.attrib['id'],
                                        precond_event.attrib.get('id'))
                        if bound_agent is not None:
                            bound_agents.append(bound_agent)
                else:
                    bound_agent = \
                        self._get_basic_agent_by_id(bound_to_term_id,
                                        precond_event.attrib.get('id'))
                    if bound_agent is not None:
                        bound_agents = [bound_agent]

            # Look for negative flag either in precondition event
            # predicate tag or in the term itself
            # (after below, neg_flag will be an object, or None)
            neg_flag = precond_event.find(
                            'predicate/mods/mod[type="ONT::NEG"]')
            negation_sign = precond_event.find('negation')
            if negation_sign is not None and negation_sign.text == '+':
                neg_flag = True
            # (after this, neg_flag will be a boolean value)
            neg_flag = neg_flag or \
                       agent_term.find('mods/mod[type="ONT::NEG"]')
            for ba in bound_agents:
                if neg_flag:
                    bc = BoundCondition(ba, False)
                else:
                    bc = BoundCondition(ba, True)
                agent.bound_conditions.append(bc)
            return
        logger.debug('Unhandled precondition event type: %s' %
                       precond_event_type)

    def _find_in_term(self, term_id, path):
        tag = self.tree.find("TERM[@id='%s']/%s" % (term_id, path))
        return tag

    def _get_basic_agent_by_id(self, term_id, event_id):
        agent = self._get_agent_by_id(term_id, event_id)
        if agent is None:
            return None
        if isinstance(agent, collections.Iterable):
            agent = agent[0]
            logger.warning('Extracting only one basic Agent from %s.'
                            % term_id)
        basic_agent = Agent(agent.name, db_refs=agent.db_refs)
        return basic_agent

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
                if residue is None:
                    residue = [None]
                if pos is None:
                    pos = [None]
                all_residues += residue
                all_pos += pos
        else:
            site_type = site_term.find("type").text
            site_name_tag = site_term.find("name")
            if site_name_tag is not None:
                site_name = site_name_tag.text
            if site_type == 'ONT::MOLECULAR-SITE':
                residue = site_term.find('features/site/code')
                if residue is not None and residue.text:
                    residue = residue.text.upper()
                else:
                    residue = None
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

    def _get_modification(self, event):
        # Find the modification type
        mod_type = event.find('type').text
        mod_type_name = mod_names.get(mod_type)
        # If it is an unknown modification type
        if mod_type_name is None:
            logger.warning('Unhandled modification type: %s')
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

    def _get_evidence(self, event_tag):
        text = self._get_evidence_text(event_tag)
        sec = self._get_section(event_tag)
        epi = {'section_type': sec}
        ev = Evidence(source_api='trips', text=text, pmid=self.doc_id,
                      epistemics=epi)
        return ev

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

    def _find_isolated_terms(self):
        all_events = self.tree.findall('EVENT')
        active_event_args = set()
        for event in all_events:
            if event.attrib.get('id') in self._static_events:
                continue
            args = event.findall('arg') + \
                   [event.find('arg1'), event.find('arg2'), event.find('arg3')]
            arg_ids = [a.attrib.get('id') for a in args if a is not None]
            active_event_args = active_event_args.union(set(arg_ids))
        all_terms = self.tree.findall('TERM')
        isolated_terms = []
        for term in all_terms:
            term_id = term.attrib.get('id')
            if term_id and term_id not in active_event_args:
                isolated_terms.append(term)
        return isolated_terms

    def _add_extracted(self, event_type, event_id):
        self.extracted_events[event_type].append(event_id)


def _get_type(element):
    type_tag = element.find('type')
    if type_tag is None:
        return None
    type_text = type_tag.text
    return type_text


def _is_type(element, type_text):
    element_type = _get_type(element)
    if element_type == type_text:
        return True
    return False


def _stmt_location_to_agents(stmt, location):
    """Apply an event location to the Agents in the corresponding Statement.

    If a Statement is in a given location we represent that by requiring all
    Agents in the Statement to be in that location.
    """
    if location is None:
        return
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


def _is_base_agent_state(agent):
    if agent.location is None and \
       not agent.mods and \
       not agent.mutations and \
       not agent.bound_conditions:
            return True
    return False


def _get_db_refs(term):
    """Extract database references for a TERM."""
    db_refs = {}
    # Here we extract the text name of the Agent
    # There are two relevant tags to consider here.
    # The <text> tag typically contains a larger phrase surrounding the
    # term but it contains the term in a raw, non-canonicalized form.
    # The <name> tag only contains the name of the entity but it is
    # canonicalized. For instance, MAP2K1 appears as MAP-2-K-1.
    agent_text_tag = term.find('name')
    if agent_text_tag is not None:
        db_refs['TEXT'] = agent_text_tag.text

    # We make a list of scored grounding terms from the DRUM terms
    grounding_terms = _get_grounding_terms(term)
    if not grounding_terms:
        # This is for backwards compatibility with EKBs without drum-term
        # scored entries. It is important to keep for Bioagents
        # compatibility.
        dbid = term.attrib.get('dbid')
        if dbid:
            dbids = dbid.split('|')
            for dbname, dbid in [d.split(':') for d in dbids]:
                if not db_refs.get(dbname):
                    db_refs[dbname] = dbid
        return db_refs, []

    # This is the INDRA prioritization of grounding name spaces. Lower score
    # takes precedence.
    ns_priority = {
        'HGNC': 1,
        'UP': 1,
        'BE': 2,
        'CHEBI': 3,
        'GO': 4,
        'FA': 5,
        'XFAM': 5,
        'NCIT': 5
    }
    # We get the top priority entry from each score group
    score_groups = itertools.groupby(grounding_terms, lambda x: x['score'])
    top_per_score_group = []
    ambiguities = []
    for score, group in score_groups:
        entries = list(group)
        for entry in entries:
            priority = 100
            for ref_ns, ref_id in entry['refs'].items():
                # Skip etc UP entries
                if ref_ns == 'UP' and ref_id == 'etc':
                    continue
                try:
                    priority = min(priority, ns_priority[ref_ns])
                except KeyError:
                    pass
                if ref_ns == 'UP':
                    if not up_client.is_human(ref_id):
                        priority = 4
            entry['priority'] = priority
        if len(entries) > 1:
            top_entry = entries[0]
            top_idx = 0
            for i, entry in enumerate(entries):
                if entry['priority'] < top_entry['priority']:
                    top_entry = entry
                    top_idx = i
            for i, entry in enumerate(entries):
                if i == top_idx:
                    continue
                if (entry['priority'] - top_entry['priority']) <= 1:
                    ambiguities.append((top_entry, entry))
        else:
            top_entry = entries[0]
        top_per_score_group.append(top_entry)
    # Get the top priority for each score group
    priorities = [entry['priority'] for entry in top_per_score_group]

    # By default, we coose the top priority entry from the highest score group
    top_grounding = top_per_score_group[0]
    # Sometimes the top grounding has much lower priority and not much higher
    # score than the second grounding. Typically 1.0 vs 0.82857 and 5 vs 2.
    # In this case we take the second entry.
    if len(top_per_score_group) > 1:
        score_diff = top_per_score_group[0]['score'] - \
                     top_per_score_group[1]['score']
        priority_diff = top_per_score_group[0]['priority'] - \
                        top_per_score_group[1]['priority']
        if score_diff < 0.2 and priority_diff >= 2:
            top_grounding = top_per_score_group[1]
    relevant_ambiguities = []
    for amb in ambiguities:
        if top_grounding not in amb:
            continue
        if top_grounding == amb[0]:
            relevant_ambiguities.append({'preferred': amb[0],
                                         'alternative': amb[1]})
        else:
            relevant_ambiguities.append({'preferred': amb[1],
                                         'alternative': amb[0]})

    for k, v in top_grounding['refs'].items():
        db_refs[k] = v

    return db_refs, relevant_ambiguities


def _get_grounding_terms(term):
    drum_terms = term.findall('drum-terms/drum-term')
    if not drum_terms:
        return None
    terms = []
    score_started = False
    for dt in drum_terms:
        # This is the primary ID
        dbid_str = dt.attrib.get('dbid')

        if not dbid_str:
            if _is_type(dt.find('types'), 'ONT::PROTEIN-FAMILY'):
                members = dt.findall('members/member')
                dbids = []
                for m in members:
                    dbid = m.attrib.get('dbid')
                    dbids.append(dbid)
                refs = {'PFAM-DEF': '|'.join(dbids)}
            # This is to handle the occasional empty drum-term
            else:
                refs = {}
        else:
            db_ns, db_id = dbid_str.split(':')
            refs = {db_ns: db_id}

        # Next look at the xref tags
        xr_tags = dt.findall('xrefs/xref')
        for xrt in xr_tags:
            dbid_str = xrt.attrib.get('dbid')
            db_ns, db_id = dbid_str.split(':')
            # XFAM xrefs are added to proteins but are
            # not desirable here
            if db_ns == 'XFAM':
                continue
            refs[db_ns] = db_id

        # Next we look at alternatives for the entry. For instance
        # we check if NCIT maps to HGNC, CHEBI, GO or BE.
        new_refs = {}
        for ref_ns, ref_id in refs.items():
            db_mappings = _get_db_mappings(ref_ns, ref_id)
            for ref_mapped in db_mappings:
                new_refs[ref_mapped[0]] = ref_mapped[1]
        for k, v in new_refs.items():
            refs[k] = v

        # Now get the match score associated with the term
        match_score = dt.attrib.get('match-score')
        db_name = dt.attrib.get('name')
        # Handling corner cases for unscored matches
        if match_score is None:
            if not score_started:
                # This is a match before other scored terms so we
                # default to 1.0
                match_score = 1.0
            else:
                # This is a match after other scored matches
                # default to a small value
                match_score = 0.1
        else:
            match_score = float(match_score)
            score_started = True
        # This is a special case to handle unscored blank drum-terms
        # at the top of the list
        if not refs:
            match_score = 0
        entity_type = dt.find('types/type')
        if entity_type is not None:
            entity_type = entity_type.text
        grounding_term = {'score': match_score,
                          'refs': refs,
                          'name': db_name,
                          'type': entity_type}
        terms.append(grounding_term)
    # Finally, the scores are sorted in descending order
    terms = sorted(terms, key=operator.itemgetter('score'), reverse=True)
    # Merge grounding terms that are identical based on the references
    # that they contain. The identical references are merged into the
    # highest scoring term. Example:
    # [{'refs': {'NCIT': '123', 'HGNC': '234'}, score: 1.0},
    # {'refs': {'HGNC': '234', 'UP', 'P123'}, score: 0.829}]
    # ==>
    # [{'refs': {'NCIT': '123', 'HGNC': '234', 'UP': 'P123'}, score: 1.0}]
    if len(terms) > 1:
        independent_terms = [terms[0]]
        for t in terms[1:]:
            any_match = False
            for it in independent_terms:
                match = False
                for k, v in t['refs'].items():
                    if k in it['refs'] and it['refs'][k] == v:
                        match = True
                        any_match = True
                if match:
                    for k, v in t['refs'].items():
                        it['refs'][k] = v
            if not any_match:
                independent_terms.append(t)
        terms = independent_terms
    return terms


def _get_db_mappings(dbname, dbid):
    db_mappings = []
    be_id = bioentities_map.get((dbname, dbid))
    if be_id is not None:
        db_mappings.append(('BE', be_id))
    if dbname == 'NCIT':
        target = ncit_map.get(dbid)
        if target is not None:
            db_mappings.append((target[0], target[1]))
    elif dbname == 'HGNC':
        standard_up_id = hgnc_client.get_uniprot_id(dbid)
        # standard_up_id will be None if the gene doesn't have a corresponding
        # protein product
        if standard_up_id:
            db_mappings.append(('UP', standard_up_id))
    elif dbname == 'UP':
        # Handle special case of UP:etc
        if not dbid == 'etc' and up_client.is_human(dbid):
            gene_name = up_client.get_gene_name(dbid)
            if gene_name:
                hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                if hgnc_id:
                    db_mappings.append(('HGNC', hgnc_id))
    return db_mappings


def _read_ncit_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../resources/ncit_map.tsv')
    ncit_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    next(csv_rows)
    for row in csv_rows:
        ncit_id = row[0]
        target_ns = row[1]
        target_id = row[2]
        ncit_map[ncit_id] = (target_ns, target_id)
    return ncit_map

ncit_map = _read_ncit_map()

def _read_bioentities_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../resources/bioentities_map.tsv')
    bioentities_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        source_ns = row[0]
        source_id = row[1]
        be_id = row[2]
        bioentities_map[(source_ns, source_id)] = be_id
    return bioentities_map

bioentities_map = _read_bioentities_map()
