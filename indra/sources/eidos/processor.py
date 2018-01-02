from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import objectpath
from indra.statements import Influence, Agent, Evidence

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

    def get_events(self):
        events = self.tree.execute("$.mentions[(@.type is 'EventMention')]")
        events = list(events)

        # Skip events that only have one argument
        #events = [e for e in events if len(e['arguments']) == 2]

        for event in events:
            # Skip events with missing arguments
            if len(event['arguments']) != 2:
                continue
            # Process causal events
            if 'Causal' in event['labels']:
                subj = event['arguments']['cause'][0]
                obj = event['arguments']['effect'][0]
            # Process origin/theme events
            if 'Origin' in event['labels']:
                subj = event['arguments']['origin'][0]
                obj = event['arguments']['theme'][0]
            subj_agent = self._get_agent(subj)
            obj_agent = self._get_agent(obj)
            subj_mods = self._get_mods(subj)
            obj_mods = self._get_mods(obj)
            # The interpretation of multiple mods is not clear yet so we
            # choose the first mod if available
            subj_delta = subj_mods[0] if subj_mods else None
            obj_delta = obj_mods[0] if obj_mods else None
            evidence = self._get_evidence(event)
            st = Influence(subj_agent, obj_agent, subj_delta, obj_delta,
                           evidence=evidence)
            self.statements.append(st)

    def _get_evidence(self, event):
        text = event.get('text')
        ev = Evidence(source_api='eidos', text=text)
        return [ev]

    def _get_mods(self, term):
        mods = []
        for mod in term.get('modifications', []):
            polarity = 'positive' if mod['type'] == 'Increase' else 'negative'
            # There is no adjective yet in the Eidos output so that is set
            # to None
            entry = {'adjective': None, 'polarity': polarity}
            mods.append(entry)
        return mods

    def _get_agent(self, term):
        name = term.get('text')
        agent = Agent(name)
        return agent
