"""
This module provides a simplified API for invoking the Turku Event Extraction
System (TEES) on text and extracting INDRA statement from TEES output.

See publication:
Jari Björne, Sofie Van Landeghem, Sampo Pyysalo, Tomoko Ohta, Filip Ginter, Yves Van de Peer, Sofia Ananiadou and Tapio Salakoski, PubMed-Scale Event Extraction for Post-Translational Modifications, Epigenetics and Protein Structural Relations. Proceedings of BioNLP 2012, pages 82-90, 2012.
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.sources.tees.processor import TEESProcessor

import os.path
import logging
import codecs

from indra.sources.tees.parse_tees import tees_parse_networkx_to_dot
import networkx.algorithms.dag as dag


logger = logging.getLogger('tees')

# If TEES isn't specified, we will check to see if any of these directories
# contain all of the files in tees_installation_files; if so, we'll assume
# that it is a TEES installation.
tees_candidate_paths = ['../TEES', '~/TEES', '~/Downloads/TEES']
tees_installation_files = ['batch.py', 'classify.py', 'train.py',
        'visualize.py']
tees_installation_dirs = ['Classifiers', 'Detectors', 'Evaluators', 'Core']

def tees_process_text(text, tees_path=None, python2_path=None):
    """Processes the specified plain text with TEES and converts output to
    supported INDRA statements.

    Parameters
    ----------
    text: str
        Plain text to process with TEES
    tees_path: str
        The path of the TEES installation directory containing classify.py.
        If None, searches several common paths.
    python2_path: str
        TEES is only compatible with python 2. This processor invokes this
        external python 2 interpreter so that the processor can be run in
        either python 2 or python 3. If None, searches for an executible named
        python2 in the PATH environment variable.

    Returns
    -------
    tp: TEESProcessor
        A TEESProcessor object which contains a list of INDRA statements
        extracted from TEES extractions
    """

    # If TEES directory is not specifies, see if any of the candidate paths
    # exist and contain all of the files expected for a TEES installation.
    for cpath in tees_candidate_paths:
        cpath = os.path.expanduser(cpath)
        #print('Checking ' + cpath)
        #print('Is it a directory?:',  os.path.isdir(cpath))
        if os.path.isdir(cpath):
            # Check to see if it has all of the expected files and directories
            has_expected_files = True
            for f in tees_installation_files:
                fpath = os.path.join(cpath, f)
                present = os.path.isfile(fpath)
                #print('Checking' , fpath , ':' , present)
                has_expected_files = has_expected_files and present

            has_expected_dirs = True
            for d in tees_installation_dirs:
                dpath = os.path.join(cpath, d)
                present = os.path.isdir(dpath)
                #print('Checking' , dpath , ':' , present)
                has_expected_dirs = has_expected_dirs and present

            if has_expected_files and has_expected_dirs:
                # We found a directory with all of the files and directories
                # we expected in a TEES installation - let's assume it's a
                # TEES installation
                tees_path = cpath
                print('Found TEES installation at' + cpath)
                break

    # If tees_path is None then we didn't find any installations
    if tees_path is None:
        raise Exception('Could not find TEES installation')

    # Try to locate python2 in one of the directories of the PATH environment
    # variable if it is not provided
    if python2_path is None:
        for path in os.environ["PATH"].split(os.pathsep):
            proposed_python2_path = os.path.join(path, 'python2.7')
            if os.path.isfile(proposed_python2_path):
                python2_path = proposed_python2_path
                print('Found python 2 interpreter at', python2_path)
                break
    if python2_path is None:
        raise Exception('Could not find python2 in the directories ' + 
                'listed in the PATH environment variable. ' + 
                'Need python2 to run TEES.')

    # Run the TEES processor
    tp = TEESProcessor(text, tees_path, python2_path)
    return tp

if __name__ == '__main__':
    #For debugging
    text = ''
    with codecs.open('abstracts_100.txt', 'r', encoding='utf-8') as f:
        text = text + f.read()
    tees_path = '/Users/daniel/Downloads/jbjorne-TEES-1125ab0'
    good_sentence = 'Raf increases the phosphorylation of BRAF.'
    weird_sentence = 'Serine 446 is constituitively phosphorylated in BRAF.'
    s = 'Our data demonstrated that Abl and Arg were activated downstream of chemokine receptors and mediated the chemokine-induced tyrosine phosphorylation of human enhancer of filamentation 1 (HEF1), an adaptor protein that is required for the activity of the guanosine triphosphatase Rap1, which mediates cell adhesion and migration.'

    tp = tees_process_text(text)
    statements = tp.statements
    for statement in statements:
        print(statement)

    # Make a graph with nodes involving Binding events
    events = tp.G
    subgraph_nodes = set()
    for node in events.node:
        if events.node[node]['is_event'] and events.node[node]['type'] == 'Gene_expression':
            subgraph_nodes.add(node)
            subgraph_nodes.update(dag.ancestors(events, node))
            subgraph_nodes.update(dag.descendants(events, node))
    print('Subgraph size: %d' % len(subgraph_nodes))
    print('Subgraph nodes: ', subgraph_nodes)
    subgraph_subset = set()
    counter = 0
    for n in subgraph_nodes:
        subgraph_subset.add(n)
        counter = counter + 1
        if counter > 200:
            break
    tees_parse_networkx_to_dot(events, 'tees_all.dot', list(events.node.keys()))

    n = tp.find_event_with_outgoing_edges('Binding', ['Theme', 'Theme2'])
    print('Events in total:', len(n))
    print('=====================')
    print('Binding events: ', n)
    print('Moo now')


    #TODO: parse binding events by looking for a Binding that has Theme and Theme2



