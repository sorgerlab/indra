from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import xml.etree.ElementTree as ET
from datetime import datetime

from indra.statements import *
from indra.statements.statements import Migration
from indra.statements.context import MovementContext
from indra.util import UnicodeXMLTreeBuilder as UTB


logger = logging.getLogger(__name__)


class CWMSError(Exception):
    pass


POLARITY_DICT = {'CC': {'ONT::CAUSE': 1,
                        'ONT::INFLUENCE': 1},
                 'EVENT': {'ONT::INCREASE': 1,
                           'ONT::MODULATE': None,
                           'ONT::DECREASE': -1,
                           'ONT::INHIBIT': -1,
                           'ONT::TRANSFORM': None,
                           'ONT::STIMULATE': 1,
                           'ONT::ARRIVE': None,
                           'ONT::DEPART': None},
                 'EPI': {'ONT::ASSOCIATE': None}}


class CWMSProcessor(object):
    """The CWMSProcessor currently extracts causal relationships between
    terms (nouns) in EKB. In the future, this processor can be extended to
    extract other types of relations, or to extract relations involving
    events.

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
    doc_id: str
        Document ID
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted from the EKB.
    sentences : dict[str: str]
        The list of all sentences in the EKB with their IDs
    paragraphs : dict[str: str]
        The list of all paragraphs in the EKB with their IDs
    par_to_sec : dict[str: str]
        A map from paragraph IDs to their associated section types
    """

    def __init__(self, xml_string):
        self.statements = []
        # Parse XML
        try:
            self.tree = ET.XML(xml_string, parser=UTB())
        except ET.ParseError:
            logger.error('Could not parse XML string')
            self.tree = None
            return

        # Get the document ID from the EKB tag.
        self.doc_id = self.tree.attrib.get('id')

        # Store all paragraphs and store all sentences in a data structure
        paragraph_tags = self.tree.findall('input/paragraphs/paragraph')
        sentence_tags = self.tree.findall('input/sentences/sentence')
        self.paragraphs = {p.attrib['id']: p.text for p in paragraph_tags}
        self.sentences = {s.attrib['id']: s.text for s in sentence_tags}
        self.par_to_sec = {p.attrib['id']: p.attrib.get('sec-type')
                           for p in paragraph_tags}

        # Keep a list of events that are part of relations and events
        # subsumed by other events
        self.relation_events = []
        self.subsumed_events = []

        # Keep a list of unhandled events for development purposes
        self._unhandled_events = []

        self._preprocess_events()

    def _preprocess_events(self):
        events = self.tree.findall("EVENT/[type]")
        for event in events:
            affected = event.find("*[@role=':AFFECTED']")
            if affected is not None:
                affected_id = affected.attrib.get('id')
                if affected_id:
                    self.subsumed_events.append(affected_id)

    def extract_causal_relations(self):
        """Extract Influence Statements from the EKB."""
        relations = self.tree.findall("CC/[type]")
        for relation in relations:
            st = self.influence_from_relation(relation)
            if st:
                self.statements.append(st)

        events = self.tree.findall("EVENT/[type]")
        for event in events:
            st = self.influence_from_event(event)
            if st:
                self.statements.append(st)
        # In some EKBs we get two redundant relations over the same arguments,
        # we eliminate these
        self._remove_multi_extraction_artifacts()

        # Print unhandled event types
        logger.debug('Unhandled event types: %s' %
                     (', '.join(sorted(list(set(self._unhandled_events))))))

    def extract_events(self):
        """Extract standalone Events from the EKB."""
        events = [(1, self.tree.findall("EVENT/[type='ONT::INCREASE']")),
                  (-1, self.tree.findall("EVENT/[type='ONT::DECREASE']"))]
        for polarity, event_list in events:
            for event_term in event_list:
                event_id = event_term.attrib.get('id')
                if event_id in self.subsumed_events or \
                        event_id in self.relation_events:
                    continue
                event = self.event_from_event(event_term)
                if event:
                    # Here we set the polarity based on the polarity implied by
                    # the increase/decrease here
                    event.delta.set_polarity(polarity)
                    self.statements.append(event)

        self._remove_multi_extraction_artifacts()

    def extract_migrations(self):
        ev_types = ['ONT::MOVE', 'ONT::DEPART', 'ONT::ARRIVE']
        events = []
        for et in ev_types:
            evs = self.tree.findall("EVENT/[type='%s']" % et)
            events += evs

        for event_term in events:
            event = self.migration_from_event(event_term)
            if event is not None:
                self.statements.append(event)

    def extract_correlations(self):
        correlations = self.tree.findall("EPI/[type='ONT::ASSOCIATE']")
        for cor in correlations:
            st = self._association_from_element(cor, 'EPI', 'NEUTRAL1',
                                                'NEUTRAL2', False)
            if st:
                self.statements.append(st)

        # self._remove_multi_extraction_artifacts()

    def _influence_from_element(self, element, element_type, subj_arg,
                                obj_arg, is_arg):
        components = self._statement_components_from_element(
            element, element_type, subj_arg, obj_arg, is_arg)
        if components is None:
            return None
        subj, obj, evidence, rel_type = components
        # If the object polarity is not given explicitly, we set it
        # based on the one implied by the relation
        if obj.delta.polarity is None:
            obj.delta.set_polarity(POLARITY_DICT[element_type][rel_type])

        st = Influence(subj, obj, evidence=[evidence])
        return st

    def influence_from_relation(self, relation):
        """Return an Influence from a CC element in the EKB."""
        return self._influence_from_element(relation, 'CC', 'FACTOR',
                                            'OUTCOME', True)

    def influence_from_event(self, event):
        """Return an Influence from an EVENT element in the EKB."""
        return self._influence_from_element(event, 'EVENT', 'AGENT',
                                            'AFFECTED', False)

    def _statement_components_from_element(self, element, element_type,
                                           member1_arg, member2_arg, is_arg):
        element_id = element.attrib.get('id')
        rel_type = element.find('type').text
        if rel_type not in POLARITY_DICT[element_type]:
            self._unhandled_events.append(rel_type)
            return None
        member1_id, member1_term = self._get_term_by_role(
            element, member1_arg, is_arg)
        member2_id, member2_term = self._get_term_by_role(
            element, member2_arg, is_arg)
        if member1_term is None or member2_term is None:
            return None

        member1 = self._get_event(member1_term)
        member2 = self._get_event(member2_term)
        if member1 is None or member2 is None:
            return None

        self.relation_events += [member1_id, member2_id, element_id]

        evidence = self._get_evidence(element)

        return member1, member2, evidence, rel_type

    def _association_from_element(self, element, element_type, member1_arg,
                                  member2_arg, is_arg):
        components = self._statement_components_from_element(
            element, element_type, member1_arg, member2_arg, is_arg)
        if components is None:
            return None
        member1, member2, evidence, _ = components
        st = Association([member1, member2], evidence=[evidence])
        return st

    def event_from_event(self, event_term):
        """Return an Event from an EVENT element in the EKB."""
        arg_id, arg_term = self._get_term_by_role(event_term, 'AFFECTED',
                                                  False)
        if arg_term is None:
            return None
        # Make an Event statement if it is a standalone event
        evidence = self._get_evidence(event_term)
        event = self._get_event(arg_term, evidence=[evidence])
        if event is None:
            return None
        event.context = self.get_context(event_term)
        return event

    def migration_from_event(self, event_term):
        """Return a Migration event from an EVENT element in the EKB."""
        # Arguments can be under AGENT or AFFECTED
        agent_arg_id, agent_arg_term = self._get_term_by_role(
            event_term, 'AGENT', False)
        affected_arg_id, affected_arg_term = self._get_term_by_role(
            event_term, 'AFFECTED', False)
        if agent_arg_term is None and affected_arg_term is None:
            return None

        # Try to get the quantitative state associated with the event
        for arg_term in [agent_arg_term, affected_arg_term]:
            if arg_term is not None:
                size_arg = arg_term.find('size')
                if size_arg:
                    break
        if size_arg is not None:
            size = self._get_size(size_arg.attrib['id'])
        else:
            size = None
        # Try to get the locations associated with the event
        neutral_id, neutral_term = self._get_term_by_role(event_term,
                                                          'NEUTRAL',
                                                          is_arg=False)
        locs = []
        if neutral_term is not None:
            text = neutral_term.find('text').text
            locs.append({'location': RefContext(name=text),
                         'role': 'origin'})

        loc = self._extract_geoloc(event_term, arg_link='to-location')
        if loc is not None:
            locs.append({'location': loc,
                         'role': 'destination'})

        loc = self._extract_geoloc(event_term, arg_link='from-location')
        if loc is not None:
            locs.append({'location': loc,
                         'role': 'origin'})

        # There are some locations at arg level, often duplicate
        if agent_arg_term:
            loc = self._extract_geoloc(agent_arg_term)
            if loc is not None:
                if loc not in [location['location'] for location in locs]:
                    locs.append({'location': loc,
                                'role': 'destination'})
        if affected_arg_term:
            loc = self._extract_geoloc(affected_arg_term)
            if loc is not None:
                if loc not in [location['location'] for location in locs]:
                    locs.append({'location': loc,
                                'role': 'destination'})

        time = self._extract_time(event_term)

        context = MovementContext(locations=locs, time=time)

        migration_grounding = 'WM/causal_factor/social_and_political/migration'
        concept = Concept('Migration',
                          db_refs={'UN': migration_grounding})
        evidence = self._get_evidence(event_term)
        event = Migration(
            concept, delta=size, context=context, evidence=[evidence])
        return event

    def _get_size(self, size_term_id):
        size_term = self.tree.find("*[@id='%s']" % size_term_id)
        value = size_term.find('value')
        if not value:
            value = size_term.find('amount')
        if value is not None:
            mod = value.attrib.get('mod')
            if mod and mod.lower() == 'almost':
                mod = 'less_than'
            value_str = value.text.strip()
            if value_str:
                value = int(value_str)
            else:
                value = None
            unit = size_term.find('unit').text.strip()
            size = QuantitativeState(value=value, unit=unit,
                                     modifier=mod)
        else:
            size = None
        return size

    def _get_term_by_role(self, term, role, is_arg):
        """Return the ID and the element corresponding to a role in a term."""
        element = term.find("%s[@role=':%s']" % ('arg/' if is_arg else '*',
                                                 role))
        if element is None:
            return None, None
        element_id = element.attrib.get('id')
        if element_id is None:
            return None, None
        element_term = self.tree.find("*[@id='%s']" % element_id)
        if element_term is None:
            return None, None
        return element_id, element_term

    def _get_event(self, event_term, evidence=None):
        """Extract and Event from the given EKB element."""
        # Now see if there is a modifier like assoc-with connected
        # to the main concept
        assoc_with = self._get_assoc_with(event_term)

        # Get the element's text and use it to construct a Concept
        element_text_element = event_term.find('text')
        if element_text_element is None:
            return None
        element_text = element_text_element.text
        element_db_refs = {'TEXT': element_text}
        element_name = sanitize_name(element_text)

        element_type_element = event_term.find('type')
        if element_type_element is not None:
            element_db_refs['CWMS'] = element_type_element.text
            # If there's an assoc-with, we tack it on as extra grounding
            if assoc_with is not None:
                element_db_refs['CWMS'] += ('|%s' % assoc_with)

        concept = Concept(element_name, db_refs=element_db_refs)

        ev_type = event_term.find('type').text
        polarity = POLARITY_DICT['EVENT'].get(ev_type)
        delta = QualitativeDelta(polarity=polarity)
        context = self.get_context(event_term)
        event_obj = Event(concept, delta=delta, context=context,
                          evidence=evidence)
        return event_obj

    def get_context(self, element):
        time = self._extract_time(element)
        geoloc = self._extract_geoloc(element)

        if time or geoloc:
            context = WorldContext(time=time, geo_location=geoloc)
        else:
            context = None
        return context

    def _extract_time(self, term):
        time = term.find('time')
        if time is None:
            time = term.find('features/time')
            if time is None:
                return None
        time_id = time.attrib.get('id')
        time_term = self.tree.find("*[@id='%s']" % time_id)
        if time_term is None:
            return None
        text = sanitize_name(time_term.findtext('text'))
        timex = time_term.find('timex')
        if timex is not None:
            year = timex.findtext('year')
            try:
                year = int(year)
            except Exception:
                year = None
            month = timex.findtext('month')
            day = timex.findtext('day')
            if year and (month or day):
                try:
                    month = int(month)
                except Exception:
                    month = 1
                try:
                    day = int(day)
                except Exception:
                    day = 1
                start = datetime(year, month, day)
                time_context = TimeContext(text=text, start=start)
            else:
                time_context = TimeContext(text=text)
        else:
            time_context = TimeContext(text=text)
        return time_context

    def _extract_geoloc(self, term, arg_link='location'):
        """Get the location from a term (CC or TERM)"""
        loc = term.find(arg_link)
        if loc is None:
            return None
        loc_id = loc.attrib.get('id')
        loc_term = self.tree.find("*[@id='%s']" % loc_id)
        if loc_term is None:
            return None
        text = loc_term.findtext('text')
        # name = loc_term.findtext('name')
        geoloc_context = RefContext(name=text)
        return geoloc_context

    def _get_assoc_with(self, element_term):
        # NOTE: there could be multiple assoc-withs here that we may
        # want to handle
        assoc_with = element_term.find('assoc-with')
        if assoc_with is not None:
            # We first identify the ID of the assoc-with argument
            assoc_with_id = assoc_with.attrib.get('id')
            # In some cases the assoc-with has no ID but has a type
            # defined in place that we can get
            if assoc_with_id is None:
                assoc_with_grounding = assoc_with.find('type').text
                return assoc_with_grounding
            # If the assoc-with has an ID then find the TERM
            # corresponding to it
            assoc_with_term = self.tree.find("*[@id='%s']" % assoc_with_id)
            if assoc_with_term is not None:
                # We then get the grounding for the term
                assoc_with_grounding = assoc_with_term.find('type').text
                return assoc_with_grounding
        return None

    def _get_evidence(self, event_tag):
        text = self._get_evidence_text(event_tag)
        sec = self._get_section(event_tag)
        epi = {'direct': False}
        if sec:
            epi['section_type'] = sec
        ev = Evidence(source_api='cwms', text=text, pmid=self.doc_id,
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

    def _remove_multi_extraction_artifacts(self):
        # Build up a dict of evidence matches keys with statement UUIDs
        evmks = {}
        logger.debug('Starting with %d Statements.' % len(self.statements))
        for stmt in self.statements:
            if isinstance(stmt, Event):
                evmk = stmt.evidence[0].matches_key()
            elif isinstance(stmt, Influence):
                evmk = (stmt.evidence[0].matches_key() +
                        stmt.subj.matches_key() + stmt.obj.matches_key())
            elif isinstance(stmt, Association):
                evmk = (stmt.evidence[0].matches_key() +
                        stmt.members[0].matches_key() +
                        stmt.members[1].matches_key())
            if evmk not in evmks:
                evmks[evmk] = [stmt.uuid]
            else:
                evmks[evmk].append(stmt.uuid)
        # This is a list of groups of statement UUIDs that are redundant
        multi_evmks = [v for k, v in evmks.items() if len(v) > 1]
        # We now figure out if anything needs to be removed
        to_remove = []
        # Remove redundant statements
        for uuids in multi_evmks:
            # Influence statements to be removed
            infl_stmts = [s for s in self.statements if (
                            s.uuid in uuids and isinstance(s, Influence))]
            infl_stmts = sorted(infl_stmts, key=lambda x: x.polarity_count(),
                                reverse=True)
            to_remove += [s.uuid for s in infl_stmts[1:]]
            # Association statements to be removed
            assn_stmts = [s for s in self.statements if (
                            s.uuid in uuids and isinstance(s, Association))]
            assn_stmts = sorted(assn_stmts, key=lambda x: x.polarity_count(),
                                reverse=True)
            # Standalone events to be removed
            events = [s for s in self.statements if (
                        s.uuid in uuids and isinstance(s, Event))]
            events = sorted(events, key=lambda x: event_delta_score(x),
                            reverse=True)
            to_remove += [e.uuid for e in events[1:]]

        # Remove all redundant statements
        if to_remove:
            logger.debug('Found %d Statements to remove' % len(to_remove))
        self.statements = [s for s in self.statements
                           if s.uuid not in to_remove]


def sanitize_name(txt):
    name = txt.replace('\n', '')
    return name


def event_delta_score(stmt):
    pol_score = 1 if stmt.delta.polarity is not None else 0
    adj_score = len(stmt.delta.adjectives)
    return (pol_score + adj_score)
