import json
from indra.statements import *
from indra.preassembler import Preassembler
from indra.ontology.world import world_ontology
from indra.preassembler.custom_preassembly import *
import indra.tools.assemble_corpus as ac


def test_event_assemble_location():
    rainfall = Concept('rainfall')
    loc1 = RefContext(name='x', db_refs={'GEOID': '1'})
    loc2 = RefContext(name='x', db_refs={'GEOID': '2'})
    ev1 = Event(rainfall, context=WorldContext(geo_location=loc1))
    ev2 = Event(rainfall, context=WorldContext(geo_location=loc2))

    pa = Preassembler(ontology=world_ontology, stmts=[ev1, ev2],
                      matches_fun=None)
    unique_stmts = pa.combine_duplicates()

    assert len(unique_stmts) == 1
    pa = Preassembler(ontology=world_ontology, stmts=[ev1, ev2],
                      matches_fun=location_matches)
    unique_stmts = pa.combine_duplicates()
    assert len(unique_stmts) == 2


def test_influence_event_hash_reference():
    rainfall = Concept('rainfall')
    loc1 = RefContext(name='x', db_refs={'GEOID': '1'})
    loc2 = RefContext(name='x', db_refs={'GEOID': '2'})
    ev1 = Event(rainfall, context=WorldContext(geo_location=loc1))
    ev2 = Event(rainfall, context=WorldContext(geo_location=loc2))
    infl = Influence(ev1, ev2)

    h1 = ev1.get_hash(refresh=True)
    h2 = ev2.get_hash(refresh=True)
    hl1 = ev1.get_hash(refresh=True, matches_fun=location_matches)
    hl2 = ev2.get_hash(refresh=True, matches_fun=location_matches)

    assert h1 == h2, (h1, h2)
    assert hl1 != hl2, (hl1, hl2)

    ij = infl.to_json(matches_fun=location_matches)
    ev1j = ev1.to_json(matches_fun=location_matches)
    assert ev1j['matches_hash'] == ij['subj']['matches_hash'],\
        (print(json.dumps(ev1j, indent=1)),
         print(json.dumps(ij, indent=1)))


def test_agent_name_custom_preassembly():
    e1 = Event(Concept('price oil'))
    e2 = Event(Concept('oil price'))
    stmts = [e1, e2]
    stmts_out = ac.run_preassembly(stmts,
                                   matches_fun=agent_name_stmt_type_matches)
    assert len(stmts_out) == 1
