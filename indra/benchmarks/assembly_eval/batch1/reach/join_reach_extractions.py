import os
import sys
import json

def join_parts(prefix):
    entities = json.load(open(prefix + '.uaz.entities.json'))
    events = json.load(open(prefix + '.uaz.events.json'))
    sentences = json.load(open(prefix + '.uaz.sentences.json'))
    full = {'events': events, 'entities': entities, 'sentences': sentences}
    return full

if __name__ == '__main__':
    prefix = sys.argv[1]
    full = join_parts(prefix)
    with open(prefix + '.json', 'wt') as fh:
        json.dump(full, fh, indent=1)
