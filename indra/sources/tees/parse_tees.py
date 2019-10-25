"""
This module invokes TEES and converts the TEES output into a networkx graph.

The top-level function is run_and_parse_tees(tees_path, text)

File format specified by:
http://2011.bionlp-st.org/home/file-formats
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible

import gzip
from lxml import etree
import tempfile
import shutil
import tarfile
import os
import os.path
import glob
import logging
import re
import subprocess
import networkx as nx
import codecs
import networkx.algorithms.dag as dag
import sys

try:  # Python 2
    basestring
except NameError:  # Python 3
    basestring = str

logger = logging.getLogger(__name__)


class TEESEntity:
    """Temporary structure for holding information about a TEES entity before
    moving this information into a networkx graph.

    Each TEESEntity corresponds to one row in the .a1 files outputted by
    TEES during extraction.

    Attributes
    ----------
    identifier : str
        The unique tag for each entity, starting with T (ex. T28)
    entity_type : str
        The type of entity (ex. Protein)
    entity_name : str
        The name of the entity, as listed in the text
    offsets : list[int]
        The lower and upper offsets where the entity was mentioned in the text
    """
    def __init__(self, identifier, entity_type, entity_name, offsets):
        self.identifier = identifier
        self.entity_type = entity_type
        self.entity_name = entity_name
        self.offsets = offsets

    def __repr__(self):
        return '%s=%s (%s, %d-%d)' % \
                (self.identifier,
                 self.entity_name,
                 self.entity_type,
                 self.offsets[0],
                 self.offsets[1])


def parse_a1(a1_text):
    """Parses an a1 file, the file TEES outputs that lists the entities in
    the extracted events.

    Parameters
    ----------
    a1_text : str
        Text of the TEES a1 output file, specifying the entities

    Returns
    -------
    entities : Dictionary mapping TEES identifiers to TEESEntity objects
        describing each entity. Each row of the .a1 file corresponds to one
        TEESEntity object.
    """
    entities = {}

    for line in a1_text.split('\n'):
        if len(line) == 0:
            continue
        tokens = line.rstrip().split('\t')
        if len(tokens) != 3:
            raise Exception('Expected three tab-seperated tokens per line ' +
                            'in the a1 file output from TEES.')

        identifier = tokens[0]
        entity_info = tokens[1]
        entity_name = tokens[2]

        info_tokens = entity_info.split()
        if len(info_tokens) != 3:
            raise Exception('Expected three space-seperated tokens in the ' + 
                            'second column of the a2 file output from TEES.')
        entity_type = info_tokens[0]
        first_offset = int(info_tokens[1])
        second_offset = int(info_tokens[2])
        offsets = (first_offset, second_offset)

        entities[identifier] = TEESEntity(
                identifier,
                entity_type,
                entity_name,
                offsets)

    return entities


def parse_a2(a2_text, entities, tees_sentences):
    """Extracts events from a TEES a2 output into a networkx directed graph.

    Parameters
    ----------
    a2_text : str
        Text of the TEES a2 file output, specifying the event graph
    sentences_xml_gz : str
        Filename with the TEES sentence segmentation in a gzipped xml format

    Returns
    -------
    events :
        A networkx graph of events. Node names are entity and event labels
        in the original A2 file (such as "E2" or "T1") and edges between nodes
        are the various properties. Text nodes (ex. "T1") have a text node
        property that gives the text.
    """
    G = nx.DiGraph()
    event_names = set()

    # Put entities into the graph
    for entity_name in entities.keys():
        offset0 = entities[entity_name].offsets[0]
        G.add_node(entity_name, text=entities[entity_name].entity_name,
                   type=entities[entity_name].entity_type, is_event=False,
                   sentence_text=tees_sentences.index_to_sentence(offset0))

    for line in a2_text.split('\n'):
        if len(line) == 0:
            continue
        if line[0] == 'T':  # New text
            tokens = line.rstrip().split('\t')
            identifier = tokens[0]
            text = tokens[2]

            if identifier not in G.nodes:
                G.add_node(identifier)
            G.nodes[identifier]['text'] = text
            G.nodes[identifier]['is_event'] = False

        elif line[0] == 'E':  # New event
            tokens = line.rstrip().split('\t')
            if len(tokens) != 2:
                raise Exception('Expected two tab-separated tokens per line ' +
                                'in TEES a2 file.')

            event_identifier = tokens[0]

            # In the second tab-separated token, we have a series of keys
            # and values separated by the colon
            key_value_pairs = tokens[1].split()
            event_name = key_value_pairs[0].split(':')[0]
            properties = dict()
            for pair in key_value_pairs:
                key_and_value = pair.split(':')
                if len(key_and_value) != 2:
                    raise Exception('Expected two colon-separated tokens ' + 
                                    'in the second column of the a2 file ' + 
                                    'output from TEES.')
                properties[key_and_value[0]] = key_and_value[1]

            # Add event to the graph if we haven't added it yet
            if event_identifier not in G.nodes:
                G.add_node(event_identifier)

            # Add edges
            for key in properties.keys():
                G.add_edge(event_identifier, properties[key],
                           relation=key)

            # We assume that node is not negated unless a event modifier
            # later says otherwise
            G.nodes[event_identifier]['negated'] = False
            G.nodes[event_identifier]['speculation'] = False
            G.nodes[event_identifier]['type'] = event_name
            G.nodes[event_identifier]['is_event'] = True

            event_names.add(event_name)

        elif line[0] == 'M':  # Event modification
            tokens = line.split('\t')
            if len(tokens) == 2:
                raise Exception('Expected two tab-separated tokens per line ' +
                                'in the a2 file output from TEES.')

            tokens2 = tokens[1].split()
            if len(tokens2) == 2:
                raise Exception('Expected two space-separated tokens per ' + 
                                'line in the a2 file output from TEES.')
            modification_type = tokens2[0]
            modified = tokens2[1]

            # But assuming this is a negation modifier, we'll need to
            # handle it
            if modification_type == 'Negation':
                G.nodes[modified]['negated'] = True
            elif modification_type == 'Speculation':
                G.nodes[modified]['speculation'] = True
            else:
                # I've only seen negation event modifiers in these outputs
                # If there are other types of modifications,
                # we'll need to handle them, since it could
                # affect whether we want to process them into statements
                print('Unknown negation event: %s' % line)
                assert(False)
    return G


class TEESSentences:
    """Parses a TEES sentences.xml.gz and creates a list of sentences and
    their corresponding indicies.

    Allows querying of a sentence corresponding to a given index.
    """

    def __init__(self, sentence_segmentations):
        self.index_to_text = dict()
        # It would be less memory intensive to use a tree, but this is simpler
        # to code

        root = etree.fromstring(sentence_segmentations.encode('utf-8'))
        for element in root.iter('sentence'):
            offset_str = element.get('charOffset')
            offset_list = offset_str.split('-')
            #
            first_offset = int(offset_list[0])
            second_offset = int(offset_list[1])
            text = element.get('text')

            for i in range(first_offset, second_offset+1):
                self.index_to_text[i] = text

    def index_to_sentence(self, index):
        """Locates a sentences in the corpus corresponding to the given index.
        If none can be found, returns None.

        Parameters
        ----------
        index : int
            Looks up a sentence with this index

        Returns
        -------
        text : str
            A sentence corresponding to the given document index, or None
            if none can be found
        """
        if index in self.index_to_text:
            return self.index_to_text[index]
        else:
            return None


def parse_output(a1_text, a2_text, sentence_segmentations):
    """Parses the output of the TEES reader and returns a networkx graph
    with the event information.

    Parameters
    ----------
    a1_text : str
        Contents of the TEES a1 output, specifying the entities
    a1_text : str
        Contents of the TEES a2 output, specifying the event graph
    sentence_segmentations : str
        Concents of the TEES sentence segmentation output XML

    Returns
    -------
    events : networkx.DiGraph
        networkx graph with the entities, events, and relationship between
        extracted by TEES
    """

    # Parse the sentence segmentation document
    tees_sentences = TEESSentences(sentence_segmentations)

    # Parse the a1 (entities) file
    entities = parse_a1(a1_text)

    # Parse the a2 (events) file
    events = parse_a2(a2_text, entities, tees_sentences)

    return events

def tees_parse_networkx_to_dot(G, output_file, subgraph_nodes):
    """Converts TEES extractions stored in a networkx graph into a graphviz
    .dot file.

    Parameters
    ----------
    G : networkx.DiGraph
        Graph with TEES extractions returned by run_and_parse_tees
    output_file : str
        Output file to which to write .dot file
    subgraph_nodes : list[str]
        Only convert the connected graph that includes these ndoes
    """

    with codecs.open(output_file, 'w', encoding='utf-8') as f:
        f.write('digraph teesParse {\n')

        mentioned_nodes = set()

        for from_node in subgraph_nodes:
            for edge in G.edges(from_node):
                to_node = edge[1]

                mentioned_nodes.add(from_node)
                mentioned_nodes.add(to_node)
                relation = G.edges[from_node, to_node]['relation']
                f.write('%s -> %s [ label = "%s" ];\n' % (from_node, to_node,
                        relation))

        for node in mentioned_nodes:
            is_event = G.nodes[node]['is_event']
            if is_event:
                node_type = G.nodes[node]['type']
                negated = G.nodes[node]['negated']
                speculation = G.nodes[node]['speculation']

                # Add a tag to the label if the event is negated or speculation
                if negated and speculation:
                    tag = ' {NS}'
                elif negated:
                    tag = ' {N}'
                elif speculation:
                    tag = ' {S}'
                else:
                    tag = ''

                node_label = node_type + tag
            else:
                node_label = G.nodes[node]['text']
            f.write('%s [label="%s"];\n' % (node, node_label))

        f.write('}\n')
