import json
from indra.statements import Activation, Agent

def process_mentions(fname):
    with open(fname, 'r') as fh:
        jd = json.load(fh)
        mentions = jd['mentions']

    events = [m for m in mentions if m['type'] == 'EventMention']
    events = [e for e in events if len(e['arguments']) == 2]
    events = [e for e in events if 'cause' in e['arguments']]

    stmts = []
    for event in events:
        cause = event['arguments']['cause'][0]['text']
        effect = event['arguments']['effect'][0]['text']
        st = Activation(Agent(cause), Agent(effect))
        stmts.append(st)
    return stmts

