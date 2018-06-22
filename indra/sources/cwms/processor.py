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

logger = logging.getLogger('cwms')


class CWMSError(Exception):
    pass


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
        self.statement_dict = {}
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
        self.extract_noun_causal_relations()
        self.extract_noun_affect_relations()
        return

    def _extract_statement(self, event):
        ev_type = _get_type(event)
        positive_ccs = {'ONT::CAUSE', 'ONT::INFLUENCE'}
        positive_events = {'ONT::INCREASE'}
        negative_events = {'ONT::DECREASE', 'ONT::INHIBIT'}
        if ev_type in positive_ccs:
            polarity = 1
            subj_arg = event.find("arg/[@role=':FACTOR']")
            obj_arg = event.find("arg/[@role=':OUTCOME']")
        elif ev_type in positive_events | negative_events:
            subj_arg = event.find("*[@role=':AGENT']")
            obj_arg = event.find("*[@role=':AFFECTED']")
            if ev_type in positive_events:
                polarity = 1
            else:
                 polarity = -1
        else:
            logger.info("Unhandled event type: %s" % ev_type)
            return Concept(ev_type)
        return self.make_statement_noun_cause_effect(event, subj_arg, obj_arg,
                                                     polarity)

    def extract_noun_causal_relations(self):
        """Extracts causal relationships between two nouns/terms (as opposed to
        events)
        """
        # Search for causal connectives of type ONT::CAUSE
        ccs = self.tree.findall("CC/[type]")
        for cc in ccs:
            self._extract_statement(cc)

    def extract_noun_affect_relations(self):
        """Extract relationships where a term/noun affects another term/noun"""
        events = self.tree.findall("EVENT/[type]")
        for event in events:
            self._extract_statement(event)

    def get_element_object(self, element):
        # Get the term with the given element id
        element_id = element.attrib.get('id')
        element_term = self.tree.find("*[@id='%s']" % element_id)

        # Not sure this is the best way to check, but it will work. It is a
        # HACKathon, afterall. Need to figure out how to propagate polarity...
        if 'EVENT' in repr(element_term):
            return self._extract_statement(element_term)

        if element_term is None:
            return

        # Get the element's text and use it to construct a Concept
        element_text_element = element_term.find('text')
        if element_text_element is None:
            return
        element_text = element_text_element.text
        element_db_refs = {'TEXT': element_text}

        element_type_element = element_term.find('type')
        if element_type_element is not None:
            element_db_refs['CWMS'] = element_type_element.text

        return Concept(element_text, db_refs=element_db_refs)

    def make_statement_noun_cause_effect(self, event_element,
                                         cause, affected, polarity):
        # Only process if both the cause and affected are present
        if cause is None or affected is None:
            return

        cause_concept = self.get_element_object(cause)
        affected_concept = self.get_element_object(affected)

        if cause_concept is None or affected_concept is None:
            return

        # Construct evidence
        ev = self._get_evidence(event_element)
        ev.epistemics['direct'] = False

        # Make statement
        if polarity == -1:
            obj_delta = {'polarity': -1, 'adjectives': []}
        else:
            obj_delta = None
        st = Influence(cause_concept, affected_concept, obj_delta=obj_delta,
                       evidence=[ev])
        if st.get_hash() not in self.statement_dict.keys():
            self.statements.append(st)
            self.statement_dict[st.get_hash()] = st
        else:
            # It's a shame we have to go this far before throwing away our work,
            # but it will do for now.
            st = self.statement_dict[st.get_hash()]
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


def _get_type(event):
    """Get the type of event."""
    children = event.getchildren()
    # What follows is another terrible hack.
    type_text_list = [ch.text for ch in children if 'type' in repr(ch)]
    if len(type_text_list) != 1:
        raise CWMSError("Unexpected event structure.")
    return type_text_list[0]
