"""This script produces a TSV that helps map the Hume ontology
to the Eidos UN ontology based on entries and examples."""

import yaml
import requests
from indra.sources import hume


def build_examples(node, tree, prefix, all_examples):
    if not prefix or prefix in ('entity', 'event', 'indicator'):
        this_prefix = node
    else:
        this_prefix = prefix + ',' + node if prefix else node
    for entry in tree:
        if isinstance(entry, str) and entry[0] != '_':
            child = entry
        elif isinstance(entry, dict):
            for child in entry.keys():
                if child[0] != '_':
                    if isinstance(entry[child], (list, dict)):
                        build_examples(child, entry[child], this_prefix,
                                       all_examples)
        else:
            continue
        if child[0] != '_':
            if this_prefix in all_examples:
                all_examples[this_prefix].add(child)
            else:
                parts = this_prefix.split(',')
                all_examples[this_prefix] = set(parts + [child])


def make_file(fname):
    hume_ont_url = ('https://raw.githubusercontent.com/BBN-E/Hume/master/'
                    'resource/ontologies/hume_ontology.yaml')

    all_examples = {}

    yml = requests.get(hume_ont_url).content
    root = yaml.load(yml)
    for top_entry in root:
        node = list(top_entry.keys())[0]
        build_examples(node, top_entry[node], None, all_examples)

    with open(fname, 'w') as fh:
        for k, v in sorted(all_examples.items(), key=lambda x: x[0]):
            fh.write('%s\t%s\n' % (k, ','.join(sorted(list(v)))))


if __name__ == '__main__':
    fname = 'hume_ontology_examples.tsv'
    make_file(fname)
