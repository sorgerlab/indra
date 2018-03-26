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
        self.statements = []
        self.extract_noun_causal_relations()

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
            # If either the factor or the outcome is missing, skip
            if factor is None or outcome is None:
                continue

            # Look up the factor term by id
            factor_id = factor.attrib.get('id')
            factor_term = self.tree.find("TERM/[@id='%s']" % factor_id)
            factor_text = factor_term.find('text').text
            factor_agent = Agent(factor_text, db_refs={'TEXT': factor_text})

            # Look up the outcome term by id
            outcome_id = outcome.attrib.get('id')
            outcome_term = self.tree.find("TERM/[@id='%s']" % outcome_id)
            outcome_text = outcome_term.find('text').text
            outcome_agent = Agent(outcome_text, db_refs={'TEXT': outcome_text})

            # Do not process if either the factor or outcome is missing
            if factor_term is None or outcome_term is None:
                continue

            # Construct evidence
            ev = self._get_evidence(cc)
            ev.epistemics['direct'] = False

            # Make statement
            st = Influence(factor_agent, outcome_agent, evidence=[ev])
            self.statements.append(st)
                
    def _get_evidence(self, event_tag):
        text = self._get_evidence_text(event_tag)
        sec = self._get_section(event_tag)
        epi = {}
        if sec:
            epi['section_type'] = sec
        ev = Evidence(source_api='cwmsreader', text=text, pmid=self.doc_id,
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

