from indra.literature import *
from indra.sources.geneways.geneways_actionmention_parser import \
        GenewaysActionMentionParser
import random
import pickle
import re
from stemming.porter2 import stem
from nltk import sent_tokenize
from nltk.tokenize import word_tokenize
from lxml import etree

from indra.resources.greek_alphabet import greek_alphabet

class FullTextMention(object):
    """Container for full text mentions and their corresponding full text"""
    def __init__(self, mention, xml_full_text):
        self.mention = mention
        self.xml_full_text = xml_full_text

    def __repr__(self):
        return '%s %s %s' % (self.mention.upstream, self.mention.actiontype,
                self.mention.downstream)

    def write_xml_to_file(self, output_file):
        text_file = open(output_file, "w")
        text_file.write(self.xml_full_text)
        text_file.close()

    def sentence_tokenize(self, text):
        #return text.split('.')
        return sent_tokenize(text)

    def get_sentences(self, root_element, block_tags):
        """Returns a list of plain-text sentences by iterating through
        XML tags except for those listed in block_tags."""
        sentences = []
        for element in root_element:
            if not self.any_ends_with(block_tags, element.tag):
                # tag not in block_tags
                if element.text is not None and not re.match('^\s*$',
                                                             element.text):
                    sentences.extend(self.sentence_tokenize(element.text))
                sentences.extend(self.get_sentences(element, block_tags))

        f = open('sentence_debug.txt', 'w')
        for s in sentences:
            f.write(s.lower() + '\n')
        f.close()
        return sentences

    def any_ends_with(self, string_list, pattern):
        """Returns true iff one of the strings in string_list ends in
        pattern."""
        try:
            s_base = basestring
        except:
            s_base = str
        is_string = isinstance(pattern, s_base)

        if not is_string:
            return False
        for s in string_list:
            if pattern.endswith(s):
                return True

        return False

    def extract_sentences(self, block_tags, strip_tags, remove_tags):
        # Remove these tags from the text
        s_text = self.xml_full_text
        for strip_tag in strip_tags:
            start_tag1 = '<' + strip_tag + '>'
            start_tag2 = '<' + strip_tag + ' [^<>]*>'
            end_tag = '</' + strip_tag + '>'
            s_text = re.sub(start_tag1, '', s_text)
            s_text = re.sub(start_tag2, '', s_text)
            s_text = re.sub(end_tag, '', s_text)

        # Remove these tags and anything in them from the text
        for remove_tag in remove_tags:
            r = '<' + remove_tag + '[^>]*>' + '[^<]*</' + remove_tag + '>'
            s_text = re.sub(r, '', s_text)

        # Convert greek characters to names
        for a in greek_alphabet.keys():
            s_text = s_text.replace(a, greek_alphabet[a])

        f = open('foo.txt', 'w')
        f.write(s_text)
        f.close()

        try:
            root = etree.fromstring(s_text.encode('utf-8'))
            sentences = self.get_sentences(root, block_tags)
        except:
            # If we failed to process xml, that probably means it's actually
            # plain text
            sentences = self.sentence_tokenize(self.xml_full_text)
        return sentences

    def find_matching_sentences(self, block_tags=None, strip_tags=None,
                                remove_tags=None):
        if block_tags is None:
            block_tags = []
        if strip_tags is None:
            strip_tags = ['italic', 'bold', 'sup', 'sub', 'xref']
        if remove_tags is None:
            remove_tags = []

        sentences = self.extract_sentences(block_tags, strip_tags, remove_tags)

        matching_sentences = []
        for sentence in sentences:
            if self.sentence_matches(sentence) or \
                self.sentence_matches(sentence.replace('-', ' ')):
                matching_sentences.append(sentence)

        return matching_sentences

    def get_tag_names(self):
        """Returns the set of tag names present in the XML."""
        root = etree.fromstring(self.xml_full_text.encode('utf-8'))
        return self.get_children_tag_names(root)

    def get_children_tag_names(self, xml_element):
        """Returns all tag names of xml element and its children."""
        tags = set()
        tags.add(self.remove_namespace_from_tag(xml_element.tag))

        for element in xml_element.iter(tag=etree.Element):
            if element != xml_element:
                new_tags = self.get_children_tag_names(element)
                if new_tags is not None:
                    tags.update(new_tags)
        return tags

    def string_matches_sans_whitespace(self, str1, str2_fuzzy_whitespace):
        """Check if two strings match, modulo their whitespace."""
        str2_fuzzy_whitespace = re.sub('\s+', '\s*', str2_fuzzy_whitespace)
        return re.search(str2_fuzzy_whitespace, str1) is not None

    def sentence_matches(self, sentence_text):
        """Returns true iff the sentence contains this mention's upstream
        and downstream participants, and if one of the stemmed verbs in
        the sentence is the same as the stemmed action type."""
        has_upstream = False
        has_downstream = False
        has_verb = False

        # Get the first word of the action type and assume this is the verb
        # (Ex. get depends for depends on)
        actiontype_words = word_tokenize(self.mention.actiontype)
        actiontype_verb_stemmed = stem(actiontype_words[0])

        words = word_tokenize(sentence_text)

        if self.string_matches_sans_whitespace(sentence_text.lower(),
            self.mention.upstream.lower()):
            has_upstream = True

        if self.string_matches_sans_whitespace(sentence_text.lower(),
            self.mention.downstream.lower()):
            has_downstream = True

        for word in words:
            if actiontype_verb_stemmed == stem(word):
                has_verb = True

        return has_upstream and has_downstream and has_verb
