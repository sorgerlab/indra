import sys
import json


def make_file(ont_json_file, fname):
    with open(ont_json_file, 'r') as fh:
        ont_json = json.load(fh)
    rows = []
    for top_key, entries in ont_json.items():
        for entry_key, examples in entries.items():
            entry_str = '%s/%s' % (top_key, entry_key)
            examples_str = '%s' % (','.join(examples))
            row = '%s\t%s' % (entry_str, examples_str)
            rows.append(row)

    with open(fname, 'w') as fh:
        fh.write('\n'.join(rows))


if __name__ == '__main__':
    ont_json_file = sys.argv[1]
    fname = 'sofia_ontology_examples.tsv'
    make_file(ont_json_file, fname)
