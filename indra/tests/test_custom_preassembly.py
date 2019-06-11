from indra.statements import *
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler.custom_preassembly import *


def test_event_assemble_location():
    rainfall = Concept('rainfall')
    loc1 = RefContext(name='x', db_refs={'GEOID': '1'})
    loc2 = RefContext(name='x', db_refs={'GEOID': '2'})
    ev1 = Event(rainfall, context=WorldContext(geo_location=loc1))
    ev2 = Event(rainfall, context=WorldContext(geo_location=loc2))

    pa = Preassembler(hierarchies=hierarchies, stmts=[ev1, ev2],
                      matches_fun=None)
    unique_stmts = pa.combine_duplicates()

    assert len(unique_stmts) == 1
    pa = Preassembler(hierarchies=hierarchies, stmts=[ev1, ev2],
                      matches_fun=location_matches)
    unique_stmts = pa.combine_duplicates()
    assert len(unique_stmts) == 2