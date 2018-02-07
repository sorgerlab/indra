from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

from indra.util import read_unicode_csv, write_unicode_csv
import logging
import argparse
import pickle
from indra.statements import Agent
import copy
import time
import math
from indra.databases.hgnc_client import get_hgnc_name

def string_is_integer(s):
    try:
        n = int(s)
        return True
    except:
        return False

class TmpDebugTestClass:
    def __init__(self, item):
        self.item = item

if __name__ == '__main__':
    logger = logging.getLogger('extract_grounding_map')

    doc = """Extracts a grounding map from a pickle of INDRA statements through
    the correspondence between the source text entity mention and the statements'
    db_refs dictionary.
    """
    doc = doc.rstrip()

    parser = argparse.ArgumentParser(description=doc)
    parser.add_argument('--input', '-i', type=str, required=True,
            help='Pickle file with a dictionary mapping each pmid to a list of ' + 
            'INDRA statements',
            dest='input_file')
    parser.add_argument('--output', '-o', type=str, required=True,
            help='Output csv test file containing the extracted grounding map',
            dest='output_file')
    args = parser.parse_args()

    # Load the statements from the pickle
    print('Loading input pickle file')
    d = pickle.load(open(args.input_file, 'rb'))
    print('\tDone!')

    # Make a list of statements - this makes it easier to track our progress
    statement_list = []
    for pmid in d:
        statement_list.extend(d[pmid])

    # Make a dictionary mapping the raw text mention to db_refs
    print('Extracting grounding information')
    text_to_refs = {}
    counter = 0
    percent_done = 0
    start_time = time.time()
    for statement in statement_list:
        # Look for attributes of the statement that are Agents
        # (The attribute name can vary from statement to statement)
        # We examine every attribute of the statement, as well as the members
        # of any attribute of the statement that is a list to search for
        # Agent objects from which to populate our grounding map
        potential_agents = []
        for a in dir(statement):
            a = getattr(statement, a)
            potential_agents.append(a)
            if hasattr(a, '__iter__'):
                for b in a:
                    potential_agents.append(b)
        for a in potential_agents:
            if isinstance(a, Agent):
                text = a.db_refs['TEXT']
                db_refs = copy.copy(a.db_refs)
                db_refs.pop('TEXT', None)

                # Convert HGNC ids to names
                if 'HGNC' in db_refs and string_is_integer(db_refs['HGNC']):
                    db_refs['HGNC'] = get_hgnc_name(db_refs['HGNC'])

                if len(db_refs.keys()) > 0:
                    text_to_refs[text] = db_refs
        counter = counter + 1

        progress = math.floor(100.0 * float(counter)/float(len(statement_list)))
        if progress > percent_done:
            percent_done = progress
            ellapsed_min = (time.time()-start_time) / 60.0
            print('%d%% done with processing statements (%f minutes ellapsed)' \
                    % (percent_done, ellapsed_min))
    print('\tDone!')

    # Convert into a list of lists
    print('Writing grounding map to file')
    refs_list = []
    for text in text_to_refs.keys():
        row = [text]
        for (db, ref) in text_to_refs[text].items():
            row.append(db)
            row.append(ref)
        refs_list.append(row)
    write_unicode_csv(args.output_file, refs_list)
    print('\tDone!')

