from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import logging
import argparse
import pickle
import copy
import time
import math
from indra.statements import Agent
from indra.databases.hgnc_client import get_hgnc_name
from indra.util import read_unicode_csv, write_unicode_csv
import indra.tools.assemble_corpus as ac


logger = logging.getLogger('extract_grounding_map')


def string_is_integer(s):
    try:
        n = int(s)
        return True
    except Exception:
        return False


if __name__ == '__main__':
    doc = """Extracts a grounding map from a pickle of INDRA statements through
    the correspondence between the source text entity mention and the
    statements' db_refs dictionary.
    """.rstrip()

    parser = argparse.ArgumentParser(description=doc)
    parser.add_argument('--input', '-i', type=str, required=True,
                        help='Pickle file with a dictionary mapping each ' +
                        'pmid to a list of INDRA statements',
                        dest='input_file')
    parser.add_argument('--output', '-o', type=str, required=True,
                        help='Output csv test file containing the extracted ' +
                        ' grounding map',
                        dest='output_file')
    args = parser.parse_args()

    # Load the statements from the pickle
    statement_list = ac.load_statements(args.input_file)

    # Make a dictionary mapping the raw text mention to db_refs
    logger.info('Extracting grounding information')
    text_to_refs = {}
    counter = 0
    percent_done = 0
    start_time = time.time()
    for statement in statement_list:
        for a in statement.agent_list():
            db_refs = copy.copy(a.db_refs)
            text = db_refs.pop('TEXT', None)

            # Convert HGNC ids to names
            if 'HGNC' in db_refs and string_is_integer(db_refs['HGNC']):
                db_refs['HGNC'] = get_hgnc_name(db_refs['HGNC'])

            if len(db_refs.keys()) > 0:
                text_to_refs[text] = db_refs
        counter = counter + 1

        progress = math.floor(100.0 * float(counter)
                              / float(len(statement_list)))
        if progress > percent_done:
            percent_done = progress
            ellapsed_min = (time.time()-start_time) / 60.0
            logger.info(('%d%% done with processing statements '
                         '(%f minutes elapsed)')
                        % (percent_done, ellapsed_min))
    logger.info('\tDone!')

    # Convert into a list of lists
    logger.info('Writing grounding map to file')
    refs_list = []
    for text in text_to_refs.keys():
        row = [text]
        for (db, ref) in text_to_refs[text].items():
            row.append(db)
            row.append(ref)
        refs_list.append(row)
    write_unicode_csv(args.output_file, refs_list)
    logger.info('\tDone!')
