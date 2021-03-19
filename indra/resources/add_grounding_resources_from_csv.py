import os
import csv
import logging
import argparse
import pandas as pd
from indra.statements.validate import validate_id

logger = logging.getLogger(__name__)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add grounding resources from'
                                     ' csv file of curations.')
    parser.add_argument('input_file')
    parser.add_argument('grounding_map_output_file')
    args = parser.parse_args()
    here = os.path.dirname(os.path.realpath(__file__))
    ignores_location = os.path.join(here, 'grounding', 'ignore.csv')
    misgrounding_location = os.path.join(here, 'grounding',
                                         'misgrounding_map.csv')
    gm_location = os.path.join(here, 'grounding',
                               args.grounding_map_output_file)
    df = pd.read_csv(args.input_file, sep=',', keep_default_na=False)
    df.loc[:, 'decision'] = df.decision.apply(lambda x: x.strip())
    logger.info('Adding ignores')
    new_ignores = df[df.decision == 'ignore']['text'].unique()
    with open(ignores_location) as f:
        old_ignores = [line.strip() for line in f.readlines()]
    with open(ignores_location, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for ignore in new_ignores:
            if ignore not in old_ignores:
                writer.writerow([ignore])
        else:
            logger.info('%s already exists in ignore.csv' % ignore)

    new_misgroundings = df[df.decision == 'misgrounding']
    new_misgroundings = new_misgroundings[['text', 'db_name', 'db_id']].\
        drop_duplicates().values

    with open(misgrounding_location, newline='') as f:
        reader = csv.reader(f, delimiter=',')
        old_rows = []
        for text, db_name, db_id in reader:
            old_rows.append((text, db_name, db_id))

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
    new_gm_entries = new_gm_entries[['text', 'decision']].drop_duplicates()
    new_gm_texts = new_gm_entries['text'].values
    new_gm_groundings = new_gm_entries['decision'].values
    new_gm_groundings = [tuple(grounding.split(':', maxsplit=1))
                         for grounding in new_gm_groundings]

    logger.info('Validating new grounding map entries by pattern')
    grounding_validity = [validate_id(db_ns, db_id)
                          for db_ns, db_id in new_gm_groundings]

    for text, (db_ns, db_id), validity in zip(new_gm_texts,
                                              new_gm_groundings,
                                              grounding_validity):
        if not validity:
            logger.warning('Invalid mapping %s,%s for text %s'
                           % (db_ns, db_id, text))

    old_rows = []
    if os.path.exists(gm_location):
        with open(gm_location, newline='') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                old_rows.append(tuple(row))
            else:
                logger.info('%s already exists in grounding map'
                            % (','.join(row)))

    with open(gm_location, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for text, (db_ns, db_id), validity in zip(new_gm_texts,
                                                  new_gm_groundings,
                                                  grounding_validity):
            if validity and (text, db_ns, db_id) not in old_rows:
                writer.writerow([text, db_ns, db_id])

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
