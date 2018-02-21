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

logger = logging.getLogger('tees_parser')


class TEESEntity:
    """Temporary structure for holding information about a TEES entity before
    moving this information into a networkx graph.

    Each TEESEntity corresponds to one row in the .a1 files outputted by
    TEES during extraction.

    Attributes
    ----------
    identifier: str
        The unique tag for each entity, starting with T (ex. T28)
    entity_type: str
        The type of entity (ex. Protein)
    entity_name: str
        The name of the entity, as listed in the text
    offsets: list[int]
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


def parse_a1(a1_filename):
    """Parses an a1 file, the file TEES outputs that lists the entities in
    the extracted events.

    Parameters
    ----------
    a1_filename: str
        File with the list of entities.

    Returns
    -------
    entities: Dictionary mapping TEES identifiers to TEESEntity objects
        describing each entity. Each row of the .a1 file corresponds to one
        TEESEntity object.
    """
    entities = {}

    with open(a1_filename, 'r') as a1_file:
        for line in a1_file:
            tokens = line.rstrip().split('\t')
            assert(len(tokens) == 3)

            identifier = tokens[0]
            entity_info = tokens[1]
            entity_name = tokens[2]

            info_tokens = entity_info.split()
            assert(len(info_tokens) == 3)
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


def parse_a2(a2_filename, entities, tees_sentences):
    """Extracts events from a TEES a2 files into a networkx directed graph.

    Parameters
    ----------
    a2_file: str
        Filename with the list of entities.
    sentences_xml_gz: str
        Filename with the TEES sentence segmentation in a gzipped xml format

    Returns
    -------
    events:
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

    with open(a2_filename, 'r') as a2_file:

        for line in a2_file:
            if line[0] == 'T':  # New text
                tokens = line.rstrip().split('\t')
                identifier = tokens[0]
                text = tokens[2]

                if identifier not in G.node:
                    G.add_node(identifier)
                G.node[identifier]['text'] = text
                G.node[identifier]['is_event'] = False

            elif line[0] == 'E':  # New event
                tokens = line.rstrip().split('\t')
                assert(len(tokens) == 2)

                event_identifier = tokens[0]

                # In the second tab-separated token, we have a series of keys
                # and values separated by the colon
                key_value_pairs = tokens[1].split()
                event_name = key_value_pairs[0].split(':')[0]
                properties = dict()
                for pair in key_value_pairs:
                    key_and_value = pair.split(':')
                    assert(len(key_and_value) == 2)
                    properties[key_and_value[0]] = key_and_value[1]

                # Add event to the graph if we haven't added it yet
                if event_identifier not in G.node:
                    G.add_node(event_identifier)

                # Add edges
                for key in properties.keys():
                    G.add_edge(event_identifier, properties[key],
                               relation=key)

                # We assume that node is not negated unless a event modifier
                # later says otherwise
                G.node[event_identifier]['negated'] = False
                G.node[event_identifier]['speculation'] = False
                G.node[event_identifier]['type'] = event_name
                G.node[event_identifier]['is_event'] = True

                event_names.add(event_name)

            elif line[0] == 'M':  # Event modification
                tokens = line.split('\t')
                assert(len(tokens) == 2)

                tokens2 = tokens[1].split()
                assert(len(tokens2) == 2)
                modification_type = tokens2[0]
                modified = tokens2[1]

                # But assuming this is a negation modifier, we'll need to
                # handle it
                if modification_type == 'Negation':
                    G.node[modified]['negated'] = True
                elif modification_type == 'Speculation':
                    G.node[modified]['speculation'] = True
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

    def __init__(self, sentences_xml_gz):
        self.index_to_text = dict()
        # It would be less memory intensive to use a tree, but this is simpler
        # to code

        with gzip.GzipFile(sentences_xml_gz, 'r') as f:
            contents = f.read()
        root = etree.fromstring(contents)
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
        index: int
            Looks up a sentence with this index

        Returns
        -------
        text: str
            A sentence corresponding to the given document index, or None
            if none can be found
        """
        if index in self.index_to_text:
            return self.index_to_text[index]
        else:
            return None


def parse_tees_output_directory(output_dir):
    """Parses the files in the TEES output directory. And returns a networkx
    graph.

    Parameters
    ----------
    output_dir: str
        Directory containing the output of the TEES system

    Returns
    -------
    events: networkx.DiGraph
        networkx graph with the entities, events, and relationship between
        extracted by TEES
    """

    # Locate the file of sentences segmented by the TEES system, described
    # in a compressed xml document
    sentences_glob = os.path.join(output_dir, '*-sentences.xml.gz')
    sentences_filename_candidates = glob.glob(sentences_glob)

    # Make sure there is exactly one such file
    if len(sentences_filename_candidates) != 1:
        m = 'Looking for exactly one file matching %s but found %d matches'
        raise Exception(m % (
            sentences_glob, len(sentences_filename_candidates)))

    # Parse the sentence segmentation document
    tees_sentences = TEESSentences(sentences_filename_candidates[0])

    # Create a temporary directory to tag the shared-task files
    tmp_dir = tempfile.mkdtemp(suffix='indra_tees_processor')

    try:
        # Make sure the tarfile with the extracted events is in shared task
        # format is in the output directory
        tarfile_glob = os.path.join(output_dir, '*-events.tar.gz')
        candidate_tarfiles = glob.glob(tarfile_glob)
        if len(candidate_tarfiles) != 1:
            raise Exception('Expected exactly one match for glob %s' %
                            tarfile_glob)

        # Decide what tar files to extract
        # (We're not blindly extracting all files because of the security
        # warning in the documentation for TarFile.extractall
        # In particular, we want to make sure that the filename doesn't
        # try to specify a relative or absolute path other than the current
        # directory by making sure the filename starts with an alphanumeric
        # character.
        # We're also only interested in files with the .a1 or .a2 extension
        tar_file = tarfile.open(candidate_tarfiles[0])
        a1_file = None
        a2_file = None
        extract_these = []
        for m in tar_file.getmembers():
            if re.match('[a-zA-Z0-9].*.a[12]', m.name):
                extract_these.append(m)

                if m.name.endswith('.a1'):
                    a1_file = m.name
                elif m.name.endswith('.a2'):
                    a2_file = m.name
                else:
                    assert(False)

        # There should be exactly two files that match these criteria
        if len(extract_these) != 2 or a1_file is None or a2_file is None:
            raise Exception('We thought there would be one .a1 and one .a2' +
                            ' file in the tarball, but we got %d files total' %
                            len(extract_these))

        # Extract the files that we decided to extract
        tar_file.extractall(path=tmp_dir, members=extract_these)

        # Parse the a1 (entities) file
        entities = parse_a1(os.path.join(tmp_dir, a1_file))

        # Parse the a2 (events) file
        events = parse_a2(os.path.join(tmp_dir, a2_file),
                          entities, tees_sentences)

        # Now that we're done, remove the temporary directory
        shutil.rmtree(tmp_dir)
        return events
    except BaseException as e:
        # If there was an exception, delete the temporary directory and
        # pass on the exception
        print('Not removing temporary directory: ' + tmp_dir)
        shutil.rmtree(tmp_dir)
        raise e


def run_and_parse_tees(text, tees_path, python2_path):
    """Runs TEES on the given text in a temporary directory and returns a
    directed networkx graph containing TEES entity and event information.

    Each node of the graph corresponds to either an entity or an event.
    Nodes have these properties:
    * is_event: True if an event, False if an entity
    * type: Specified only if the node is an event; gives the even type (ex.
        "Phosphorylation")
    * text: Specified only if the node is an entity; gives the text describing
        the entity in the original plain text (ex. "BRAF"), rather than using
        some standardized identifier

    Edges have the property relation, that describe the relationship between
    two nodes as listed in the originl .a2 file.

    Invokes TEES by calling a new python interpreter so that although TEES
    is only compatable with python 2, this script can be used with either
    python 2 or python 3.

    Parameters
    ----------
    text: str
        Text from which to extract relationships
    tees_path: str
        Path to the TEES directory
    python2_path: str
        The path to the python 2 interpreter

    Returns
    -------
    events: networkx.DiGraph
        networkx graph containing the entities, events, and relationships
        extracted by TEES
    """
    # Make sure the provided TEES directory exists
    if not os.path.isdir(tees_path):
        raise Exception('Provided TEES directory does not exist.')

    # Make sure the classify.py script exists within this directory
    classify_path = 'classify.py'
    # if not os.path.isfile(classify_path):
    #    raise Exception('classify.py does not exist in provided TEES path.')

    # Create a temporary directory to tag the shared-task files
    tmp_dir = tempfile.mkdtemp(suffix='indra_tees_processor')

    pwd = os.path.abspath(os.getcwd())

    try:
        # Write text to a file in the temporary directory
        text_path = os.path.join(tmp_dir, 'text.txt')
        # Had some trouble with non-ascii characters. A possible TODO item in
        # the future is to look into resolving this, for now just ignoring
        # non-ascii characters
        with codecs.open(text_path, 'w', encoding='ascii', errors='ignore') \
                as f:
            f.write(text)

        # Run TEES
        output_path = os.path.join(tmp_dir, 'output')
        model_path = os.path.join(tees_path, 'tees_data/models/GE11-test/')
        command = [python2_path, classify_path, '-m', model_path,
                   '-i', text_path,
                   '-o', output_path]
        try:
            pwd = os.path.abspath(os.getcwd())
            os.chdir(tees_path)  # Change to TEES directory
            # print('cwd is:', os.getcwd())
            # out = subprocess.check_output(command, stderr=subprocess.STDOUT)
            p = subprocess.Popen(command, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, cwd=tees_path)
            p.wait()
            (so, se) = p.communicate()
            print(so)
            print(se)
            os.chdir(pwd)  # Change back to previous directory
            # print('cwd is:', os.getcwd())
            # print(out.decode('utf-8'))
        except BaseException as e:
            # If there's an error, print it out and then propagate the
            # exception
            os.chdir(pwd)  # Change back to previous directory
            # print (e.output.decode('utf-8'))
            raise e

        # Parse TEES output
        events = parse_tees_output_directory(tmp_dir)

        # Remove the temorary directory
        shutil.rmtree(tmp_dir)

    except BaseException as e:
        # If there was an exception, delete the temporary directory and
        # pass on the exception
        shutil.rmtree(tmp_dir)
        raise e
    # Return parsed TEES output
    # print('Events: ' , events)
    return events


def tees_parse_networkx_to_dot(G, output_file, subgraph_nodes):
    """Converts TEES extractions stored in a networkx graph into a graphviz
    .dot file.

    Parameters
    ----------
    G: networkx.DiGraph
        Graph with TEES extractions returned by run_and_parse_tees
    output_file: str
        Output file to which to write .dot file
    subgraph_nodes: list[str]
        Only convert the connected graph that includes these ndoes
    """

    with codecs.open(output_file, 'w', encoding='utf-8') as f:
        f.write('digraph teesParse {\n')

        mentioned_nodes = set()

        for from_node in subgraph_nodes:
            for to_node in G.edge[from_node]:
                mentioned_nodes.add(from_node)
                mentioned_nodes.add(to_node)
                relation = G.edge[from_node][to_node]['relation']
                f.write('%s -> %s [ label = "%s" ];\n' % (from_node, to_node,
                        relation))

        for node in mentioned_nodes:
            is_event = G.node[node]['is_event']
            if is_event:
                node_type = G.node[node]['type']
                negated = G.node[node]['negated']
                speculation = G.node[node]['speculation']

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
                node_label = G.node[node]['text']
            f.write('%s [label="%s"];\n' % (node, node_label))

        f.write('}\n')
