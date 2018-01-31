import sys
import json
import time
import pickle
from indra.sources.trips import process_xml
from indra.sources.trips.drum_reader import DrumReader


def set_pmid(statements, pmid):
    for stmt in statements:
        for evidence in stmt.evidence:
            evidence.pmid = pmid


def read_content(content, host):
    all_statements = []
    for pmid, sentences in content.items():
        print('================================')
        print('Processing %d sentences for %s' % (len(sentences), pmid))
        ts = time.time()
        dr = DrumReader(to_read=sentences, host=host)
        try:
            dr.start()
        except SystemExit:
            pass
        statements = []
        for extraction in dr.extractions:
            tp = process_xml(extraction)
            statements += tp.statements
        set_pmid(statements, pmid)
        te = time.time()
        print('Reading took %d seconds and produced %d Statements.' %
              (te-ts, len(statements)))
        all_statements += statements
    return all_statements


def save_results(statements, out_fname):
    with open(out_fname, 'wb') as fh:
        pickle.dump(statements, fh)


if __name__ == '__main__':
    host = sys.argv[1]
    file_name = sys.argv[2]
    with open(file_name, 'rt') as fh:
        content = json.load(fh)
    statements = read_content(content, host)
    save_results(statements, 'results.pkl')
