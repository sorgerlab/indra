from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import xml.etree.ElementTree as ET
from indra.statements import *
from indra.util import UnicodeXMLTreeBuilder as UTB

logger = logging.getLogger('cwms')


class CWMSError(Exception):
    pass


POLARITY_DICT = {'CC': {'ONT::CAUSE': 1,
                        'ONT::INFLUENCE': 1},
                 'EVENT': {'ONT::INCREASE': 1,
                           'ONT::MODULATE': None,
                           'ONT::DECREASE': -1,
                           'ONT::INHIBIT': -1,
                           'ONT::TRANSFORM': None,
                           'ONT::STIMULATE': 1}}


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

        # Extract statements
        self.extract_noun_relations('CC')
        self.extract_noun_relations('EVENT')
        return

    def _get_subj_obj(self, event):
        """Get the concepts for a relation given and element.

        The ontological type of the event is used to infer the labels of agents
        and the polarity of the influence (see `_positive_ccs`,
        `_positive_events`, and `_negative_events` class attributes).
        """
        ev_type = event.find('type').text
        if ev_type in POLARITY_DICT['CC'].keys():
            polarity = POLARITY_DICT['CC'][ev_type]
            subj = self._get_concept(event, "arg/[@role=':FACTOR']")
            obj = self._get_concept(event, "arg/[@role=':OUTCOME']")
        elif ev_type in POLARITY_DICT['EVENT'].keys():
            polarity = POLARITY_DICT['EVENT'][ev_type]
            subj = self._get_concept(event, "*[@role=':AGENT']")
            obj = self._get_concept(event, "*[@role=':AFFECTED']")
        else:
            logger.debug("Unhandled event type: %s" % ev_type)
            return None, None, None

        return subj, obj, polarity

    def extract_noun_relations(self, key):
        """Extract relationships where a term/noun affects another term/noun"""
        events = self.tree.findall("%s/[type]" % key)
        for event in events:
            subj, obj, pol = self._get_subj_obj(event)
            self._make_statement_noun_cause_effect(event, subj, obj, pol)

    def _get_concept(self, event, find_str):
        """Get a concept referred from the event by the given string."""
        # Get the term with the given element id
        element = event.find(find_str)
        if element is None:
            return
        element_id = element.attrib.get('id')
        element_term = self.tree.find("*[@id='%s']" % element_id)

        if element_term is None:
            return

        # Now see if there is a modifier like assoc-with connected
        # to the main concept
        assoc_with = self._get_assoc_with(element_term)

        # Get the element's text and use it to construct a Concept
        element_text_element = element_term.find('text')
        if element_text_element is None:
            return
        element_text = element_text_element.text
        element_db_refs = {'TEXT': element_text}

        element_type_element = element_term.find('type')
        if element_type_element is not None:
            element_db_refs['CWMS'] = element_type_element.text
            # If there's an assoc-with, we tack it on as extra grounding
            if assoc_with is not None:
                element_db_refs['CWMS'] += ('|%s' % assoc_with)

        return Concept(element_text, db_refs=element_db_refs)

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

    def _make_statement_noun_cause_effect(self, event_element,
                                          cause_concept, affected_concept,
                                          polarity):
        """Make the Influence statement from the component parts."""
        if cause_concept is None or affected_concept is None:
            return

        # Construct evidence
        ev = self._get_evidence(event_element)
        ev.epistemics['direct'] = False

        # Make statement
        obj_delta = {'polarity': polarity, 'adjectives': []}
        st = Influence(cause_concept, affected_concept, obj_delta=obj_delta,
                       evidence=[ev])
        self.statements.append(st)
        return st

    def _get_evidence(self, event_tag):
        text = self._get_evidence_text(event_tag)
        sec = self._get_section(event_tag)
        epi = {}
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
