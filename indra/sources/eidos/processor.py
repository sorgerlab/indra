from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import objectpath
from indra.statements import Activation, Agent

class EidosProcessor(object):
    """The EidosProcessor extracts INDRA Statements from Eidos output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Eidos extractions.

    Attributes
    ----------
    tree : objectpath.Tree
        The objectpath Tree object representing the extractions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    all_events : dict[str, str]
        The frame IDs of all events by type in the Eidos extraction.
    """
    def __init__(self, json_dict):
        self.tree = objectpath.Tree(json_dict)
        self.statements = []

    def get_events(file_name):
        with open(file_name, 'r') as fh:
            jd = json.load(fh)
            mentions = jd['mentions']

        # Just extract event mentions
        events = [m for m in mentions if m['type'] == 'EventMention']
        # Skip events that only have one argument
        events = [e for e in events if len(e['arguments']) == 2]

        stmts = []
        for event in events:
            if 'Causal' in event['labels']:
                cause = event['arguments']['cause'][0]['text']
                effect = event['arguments']['effect'][0]['text']
                st = Activation(Agent(cause), Agent(effect))
                stmts.append(st)
            if 'Origin' in event['labels']:
                origin = event['arguments']['origin'][0]['text']
                theme = event['arguments']['theme'][0]['text']
                st = Activation(Agent(origin), Agent(theme))
                stmts.append(st)
        return stmts

