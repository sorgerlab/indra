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
        self.extract_noun_causal_relations()
        self.extract_noun_inhibit_relations()
        self.extract_noun_influence_relations()

    def extract_noun_causal_relations(self):
        """Extracts causal relationships between two nouns/terms (as opposed to
        events)
        """
        # Search for causal connectives of type ONT::CAUSE
        ccs = self.tree.findall("CC/[type='ONT::CAUSE']")
        for cc in ccs:
            # Each cause should involve a factor term and an outcome term
            factor = cc.find("arg/[@role=':FACTOR']")
            outcome = cc.find("arg/[@role=':OUTCOME']")
            polarity = 1
            self.make_statement_noun_cause_effect(cc, factor, outcome,
                                                  polarity)

    def extract_noun_inhibit_relations(self):
        """Extracts relationships involving one term/noun inhibiting another
        term/noun"""
        events = self.tree.findall("EVENT/[type='ONT::INHIBIT']")
        events += self.tree.findall("EVENT/[type='ONT::DECREASE']")
        for event in events:
            # Each inhibit event should involve an agent and an affected
            agent = event.find("*[@role=':AGENT']")
            affected = event.find("*[@role=':AFFECTED']")
            polarity = -1
            self.make_statement_noun_cause_effect(event, agent, affected,
                                                  polarity)

    def extract_noun_influence_relations(self):
        """Extracts relationships with one term influencing another term."""
        ccs = self.tree.findall("CC/[type='ONT::INFLUENCE']")
        for cc in ccs:
            # Each cause should involve a factor term and an outcome term
            factor = cc.find("arg/[@role=':FACTOR']")
            outcome = cc.find("arg/[@role=':OUTCOME']")
            polarity = 1
            self.make_statement_noun_cause_effect(cc, factor, outcome,
                                                  polarity)

    def make_statement_noun_cause_effect(self, event_element,
                                         cause, affected, polarity):
        # Only process if both the cause and affected are present
        if cause is None or affected is None:
            return

        # Get the term with the given cause id
        cause_id = cause.attrib.get('id')
        cause_term = self.tree.find("TERM/[@id='%s']" % cause_id)
        if cause_term is None:
            return

        # Get the cause's text and use it to construct a Concept
        cause_text_element = cause_term.find('text')
        if cause_text_element is None:
            return
        cause_text = cause_text_element.text
        cause_db_refs = {'TEXT': cause_text}
        #
        cause_type_element = cause_term.find('type')
        if cause_type_element is not None:
            cause_db_refs['CWMS'] = cause_type_element.text
        #
        cause_concept = Concept(cause_text, db_refs=cause_db_refs)

        # Get the term with the given affected id
        affected_id = affected.attrib.get('id')
        affected_term = self.tree.find("TERM/[@id='%s']" % affected_id)
        if affected_term is None:
            return

        # Get the affected's text and type and use them to construct a Concept
        affected_text_element = affected_term.find('text')
        if affected_text_element is None:
            return
        affected_text = affected_text_element.text
        affected_db_refs = {'TEXT': affected_text}
        #
        affected_type_element = affected_term.find('type')
        if affected_type_element is not None:
            affected_db_refs['CWMS'] = affected_type_element.text
        #
        affected_concept = Concept(affected_text, db_refs=affected_db_refs)

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
        self.statements.append(st)

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
