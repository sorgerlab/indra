from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
import logging
import objectpath
from indra.statements import Influence, Agent, Evidence


logger = logging.getLogger('eidos')


class EidosJsonLdProcessor(object):
    """This processor extracts INDRA Statements from Eidos JSON-LD output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Eidos extractions in JSON-LD format.

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
        events = \
            self.tree.execute("$.extractions[(@.@type is 'DirectedRelation')]")
        if not events:
            return

        entities = \
            self.tree.execute("$.extractions[(@.@type is 'Entity')]")
        entity_ids = \
            self.tree.execute("$.extractions[(@.@type is 'Entity')].@id")
        entity_dict = {id:entity for id, entity in zip(entity_ids, entities)}

        # The first state corresponds to increase/decrease
        def get_polarity(x):
            # x is either subj or obj
            if 'states' in x.keys():
                if x['states'][0]['type'] == 'DEC':
                    return -1
                elif x['states'][0]['type'] == 'INC':
                    return 1
                else:
                    return None
            else:
                return None

        def get_adjectives(x):
            # x is either subj or obj
            if 'states' in x.keys():
                if 'modifiers' in x['states'][0].keys():
                    return [mod['text'] for mod in
                            x['states'][0]['modifiers']]
            else:
                return []

        def _get_eidos_groundings(entity):
            """Return Eidos groundings are a list of tuples with
            scores."""
            return [(g['ontologyConcept'], g['value'])
                    for g in entity.get('grounding', [])]

        def _make_agent(entity):
            """Return an Agent from an Eidos entity."""
            # For now we just use the text for the agent as the name
            name = entity['text']
            # Save raw text and Eidos scored groundings as db_refs
            db_refs = {'TEXT': entity['text'],
                    'EIDOS': _get_eidos_groundings(entity)}
            agent = Agent(name, db_refs=db_refs)
            return agent

        for event in events:
            if 'Causal' in event['labels']:
                # For now, just take the first source and first destination.
                # Later, might deal with hypergraph representation.
                subj = entity_dict[event['sources'][0]['@id']]
                obj = entity_dict[event['destinations'][0]['@id']]

                subj_delta = {'adjectives': get_adjectives(subj),
                            'polarity': get_polarity(subj)}
                obj_delta = {'adjectives': get_adjectives(obj),
                            'polarity': get_polarity(obj)}

                evidence = self._get_evidence(event)

                st = Influence(_make_agent(subj), _make_agent(obj),
                            subj_delta, obj_delta, evidence=evidence)

                self.statements.append(st)

    @staticmethod
    def _get_evidence(event):
        text = event.get('text')
        annotations = {
                'found_by'   : event.get('rule'),
                'provenance' : event.get('provenance'),
                }
        ev = Evidence(source_api='eidos', text=text, annotations=annotations)
        return [ev]


class EidosJsonProcessor(object):
    """This processor extracts INDRA Statements from Eidos JSON (not JSON-LD)
    output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the Eidos extractions in JSON
        (not JSON-LD) format.

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
            elif 'Origin' in event['labels']:
                subj = event['arguments']['origin'][0]
                obj = event['arguments']['theme'][0]
            # Skip correlation events for now
            elif 'Correlation' in event['labels']:
                logger.warning('Correlation event %s skipped.' % event['id'])
                continue
            else:
                logger.warning('Could not classify event with labels: %s' %
                               ', '.join(event['labels']))
                continue
            subj_agent = self._get_agent(subj)
            obj_agent = self._get_agent(obj)
            subj_mods = self._get_mods(subj)
            obj_mods = self._get_mods(obj)
            # The interpretation of multiple mods is not clear yet so we
            # choose the first mod if available
            subj_delta = subj_mods[0] if subj_mods else \
                {'adjectives': [], 'polarity': None}
            obj_delta = obj_mods[0] if obj_mods else \
                {'adjectives': [], 'polarity': None}
            evidence = self._get_evidence(event)
            st = Influence(subj_agent, obj_agent, subj_delta, obj_delta,
                           evidence=evidence)
            self.statements.append(st)

    @staticmethod
    def _get_evidence(event):
        text = event.get('text')
        annotations = {'found_by' : event['foundBy']}
        ev = Evidence(source_api='eidos', text=text, annotations=annotations)
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
