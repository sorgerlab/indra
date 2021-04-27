import os
import json
from indra.sources import eidos
from indra.sources.eidos.bio_processor import get_agent_bio
from indra.statements import Concept


path_this = os.path.dirname(os.path.abspath(__file__))


def test_process_text_bio():
    ep = eidos.process_text_bio('virus increases death')
    assert ep is not None
    assert len(ep.statements) == 1
    stmt = ep.statements[0]
    from indra.statements import Activation
    assert isinstance(stmt, Activation)


def test_sanitize():
    # Make sure sanitization works
    sanitized = eidos.processor._sanitize('-LRB-something-RRB-')
    assert sanitized == '(something)'


def test_bio_entity_extract():
    jsonld = os.path.join(path_this, 'eidos_bio_abstract.json')
    with open(jsonld, 'r') as fh:
        js = json.load(fh)
    agents = eidos.process_json_bio_entities(js)
    assert len(agents) == 11
    from indra.statements import Agent
    assert all(isinstance(a, Agent) for a in agents)
    ag = [a for a in agents if a.name == 'Therapeutics'][0]
    assert ag.db_refs['MESH'] == 'D013812'
    assert ag.db_refs['EFO'] == '0000727'


def test_get_agent_bio():
    # (raw text, normalized text, groundings, name)
    groundings = (
        ('xxx', 'yyy', {}, 'yyy'),
        ('xxx', 'checklist', {'MESH': 'D057189'}, 'Checklist'),
        ('checklist', 'yyy', {'MESH': 'D057189'}, 'Checklist'),
        ('checklist', 'life insurance', {'MESH': 'D057189'}, 'Checklist')
    )

    for raw_text, norm_text, groundings, name in groundings:
        concept = Concept(norm_text, db_refs={'TEXT': raw_text})
        agent = get_agent_bio(concept)
        assert agent.name == name, agent
        for ns, id in groundings.items():
            assert agent.db_refs.get(ns) == id, agent.db_refs
        assert agent.db_refs['TEXT'] == raw_text
        assert agent.db_refs['TEXT_NORM'] == norm_text


def test_bio_custom_grounding():
    def my_grounder(txt, context):
        return {'MYDB': 'MYGROUNDING'}
    agent = get_agent_bio(Concept('x',
                                  {'TEXT': 'x'}),
                          grounder=my_grounder)
    assert agent.db_refs['MYDB'] == 'MYGROUNDING'
