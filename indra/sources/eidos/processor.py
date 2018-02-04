from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import logging
import objectpath
from indra.statements import Influence, Agent, Evidence


logger = logging.getLogger('eidos')


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
            subj_delta = subj_mods[0] if subj_mods else {'adjectives': None, 'polarity': None}
            obj_delta = obj_mods[0] if obj_mods else {'adjectives':None, 'polarity': None}
            evidence = self._get_evidence(event)
            st = Influence(subj_agent, obj_agent, subj_delta, obj_delta,
                           evidence=evidence)
            self.statements.append(st)

    @staticmethod
    def _get_evidence(event):
        text = event.get('text')
        ev = Evidence(source_api='eidos', text=text)
        return [ev]

    @staticmethod
    def _get_mods(term):
        mods = []
        attachments = term.get('attachments', [])
        if len(attachments) > 1:
            logger.warning('More than one attachment to event.')
        for attachment in attachments:
            # Get the polarity
            if attachment['type'] == 'Increase':
                polarity = 1
            elif attachment['type'] == 'Decrease':
                polarity = -1
            else:
                polarity = None
            # Get the adjective
            mod = attachment.get('mod')
            mod_dict = json.loads(mod)
            adjectives = mod_dict.get('quantifier', [])
            entry = {'adjectives': adjectives, 'polarity': polarity}
            mods.append(entry)
        return mods

    @staticmethod
    def _get_agent(term):
        name = term.get('text')
        agent = Agent(name)
        return agent
