"""This script allows users to add grounding resources from a csv file.

The script can add entries to ignore.csv, misgrounding_map.csv, and grounding
map entries either to an existing grounding map file in
indra/resources/grounding/ or to a new grounding map file. The script will not
add entries that already exist in a resource file. New grounding map entries
are checked to ensure they have a valid namespace prefix and identifier. See
the help text for this script for info on the proper format for input csv files.
"""

import os
import csv
import logging
import argparse
import pandas as pd
from indra.statements.validate import validate_id

logger = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add grounding resources from'
                                     ' csv file of curations. Can add entries'
                                     ' to ignore.csv, misgrounding_map.csv'
                                     ' and grounding_map entries to new or'
                                     ' existing grounding map file in'
                                     ' indra/resources/grounding/. Does not'
                                     ' add entries that already exist in these'
                                     ' files. Checks that new grounding map'
                                     ' entries have valid namespace prefix,'
                                     ' identifier pairs.')

    parser.add_argument('input_file',
                        help='CSV file containing curated groundings. It must'
                        ' have columns "text", "db_name", "db_id", and'
                        ' "decision". For a given row, the "text" entry'
                        ' contains a raw agent text. "db_name" and "db_id"'
                        ' contain a namespace prefix and ID to which the agent'
                        ' text has been mapped by some reading system. The'
                        ' "decision" column specifies what should done,'
                        ' if anything, to correct grounding for this row.'
                        ' Valid values in the decision column are "ignore"'
                        ' for adding the agent text for the row to'
                        ' ignore.csv, "misgrounding" for adding the'
                        ' original grounding db_name,db_id as a misgrounding'
                        ' for the agent text, f"{db_name2}:{db_id2}" to map'
                        ' the agent text to db_name2,db_id2 in the grounding'
                        ' map. Other valid entries are "correct" and "other",'
                        ' to specify that the original grounding is correct'
                        ' or that some other action is needed beyond making'
                        ' an entry in the grounding resource files. Warnings'
                        ' will be emitted for rows with invalid decisions.')

    parser.add_argument('grounding_map_output_file',
                        help='Base filename for placing new grounding'
                        ' map entries. File will be located in directory'
                        ' "indra/resources/grounding/". A new file will'
                        ' be created if specified filename does not exist'
                        ' in that directory. Otherwise entries will be'
                        ' appended to the existing file.')

    args = parser.parse_args()
    # get paths to resource files
    here = os.path.dirname(os.path.realpath(__file__))
    ignores_location = os.path.join(here, 'grounding', 'ignore.csv')
    misgrounding_location = os.path.join(here, 'grounding',
                                         'misgrounding_map.csv')
    gm_location = os.path.join(here, 'grounding',
                               args.grounding_map_output_file)
    df = pd.read_csv(args.input_file, sep=',', keep_default_na=False,
                     usecols=['text', 'db_name', 'db_id', 'decision'],
                     dtype=str)
    # strip extra whitespace from decision column to account for
    # common curation error
    df.loc[:, 'decision'] = df.decision.apply(lambda x: x.strip())
    logger.info('Adding ignores')
    # Filter to unique values to avoid adding entries multiple times
    new_ignores = df[df.decision == 'ignore']['text'].unique()
    # Gather existing ignores so we can avoid adding any of them again
    with open(ignores_location) as f:
        old_ignores = [line.strip() for line in f.readlines()]
    # Add new entries
    with open(ignores_location, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for ignore in new_ignores:
            if ignore not in old_ignores:
                writer.writerow([ignore])
            else:
                logger.info('%s already exists in ignore.csv' % ignore)
    
    new_misgroundings = df[df.decision == 'misgrounding']
    # Filter duplicates to avoid adding same entry multiple times
    new_misgroundings = new_misgroundings[['text', 'db_name', 'db_id']].\
        drop_duplicates().values

    # Gather existing entries to avoid adding again
    with open(misgrounding_location, newline='') as f:
        reader = csv.reader(f, delimiter=',')
        old_rows = []
        for text, db_name, db_id in reader:
            old_rows.append((text, db_name, db_id))
    # Add new entries
    with open(misgrounding_location, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for text, db_name, db_id in new_misgroundings:
            if (text, db_name, db_id) not in old_rows:
                writer.writerow([text, db_name, db_id])
            else:
                logger.info('%s,%s,%s already exists in misgrounding_map.csv'
                            % (text, db_name, db_id))
    # Filter out Famplex entries since they belong in the Famplex grounding map
    new_gm_entries = df[df.decision.apply(lambda x: ':' in x
                                          and not x.startswith('FPLX'))]
    # Remove duplicate new entries
    new_gm_entries = new_gm_entries[['text', 'decision']].drop_duplicates()
    new_gm_texts = new_gm_entries['text'].values
    new_gm_groundings = new_gm_entries['decision'].values
    new_gm_groundings = [tuple(grounding.split(':', maxsplit=1))
                         for grounding in new_gm_groundings]
    # Check validity of new grounding map entries
    logger.info('Validating new grounding map entries by pattern')
    grounding_validity = [validate_id(db_ns, db_id)
                          for db_ns, db_id in new_gm_groundings]

    for text, (db_ns, db_id), validity in zip(new_gm_texts,
                                              new_gm_groundings,
                                              grounding_validity):
        if not validity:
            logger.warning('Invalid mapping %s,%s for text %s'
                           % (db_ns, db_id, text))
    # Gather old entries to avoid adding again
    old_rows = []
    if os.path.exists(gm_location):
        with open(gm_location, newline='') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                old_rows.append(tuple(row))

    # Add new entries
    with open(gm_location, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for text, (db_ns, db_id), validity in zip(new_gm_texts,
                                                  new_gm_groundings,
                                                  grounding_validity):
            if (text, db_ns, db_id) not in old_rows:
                writer.writerow([text, db_ns, db_id])
            elif validity:
                logger.info('%s,%s,%s already exists in grounding map'
                            % (text, db_ns, db_id))
    # Gather invalid rows and emit warnings for them
    invalid_rows = df[~((df.decision == 'other') |
                        (df.decision == 'correct') |
                        (df.decision == 'ignore') |
                        (df.decision == 'misgrounding') |
                        (df.decision.apply(lambda x: ':' in x)) |
                        (df.decision.apply(lambda x: not x)))]

    for _, row in invalid_rows.iterrows():
        logger.warning('Invalid decision for '
                       '%s,%s,%s: %s' % (row.text, row.db_name,
                                         row.db_id, row.decision))
