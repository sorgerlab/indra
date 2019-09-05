from indra.statements import *


def has_location(stmt):
    """Return True if a Statement has grounded geo-location context."""
    if isinstance(stmt, Migration):
        if not stmt.context or not stmt.context.locations:
            return False
    elif not stmt.context or not stmt.context.geo_location or \
        not (stmt.context.geo_location.db_refs.get('GEOID') or
             stmt.context.geo_location.name):
        return False
    return True


def has_time(stmt):
    """Return True if a Statement has time context."""
    if not stmt.context or not stmt.context.time or \
            stmt.context.time.text == "None":
        return False
    return True


def get_location_from_object(loc_obj):
    """Return geo-location from a RefContext location object."""
    if loc_obj.db_refs.get('GEOID'):
        return loc_obj.db_refs['GEOID']
    elif loc_obj.name:
        return loc_obj.name
    else:
        return None


def get_location(stmt):
    """Return the grounded geo-location context associated with a Statement."""
    if not has_location(stmt):
        location = None
    elif isinstance(stmt, Migration):
        location = []
        for loc in stmt.context.locations:
            loc_obj = loc['location']
            location.append((get_location_from_object(loc_obj), loc['role']))
    else:
        location = get_location_from_object(stmt.context.geo_location)
    return location


def get_time(stmt):
    """Return the time context associated with a Statement."""
    if not has_time(stmt):
        time = None
    else:
        time = stmt.context.time
    return time


def location_matches(stmt):
    """Return a matches_key which takes geo-location into account."""
    if isinstance(stmt, Event):
        context_key = get_location(stmt)
        matches_key = str((stmt.concept.matches_key(), context_key))
    elif isinstance(stmt, Influence):
        subj_context_key = get_location(stmt.subj)
        obj_context_key = get_location(stmt.obj)
        matches_key = str((stmt.matches_key(), subj_context_key,
                           obj_context_key))
    else:
        matches_key = stmt.matches_key()
    return matches_key


def event_location_refinement(st1, st2, hierarchies):
    """Return True if there is a location-aware refinement between Events."""
    ref = st1.refinement_of(st2, hierarchies)
    if not ref:
        return False
    if not has_location(st2):
        return True
    elif not has_location(st1):
        return False
    else:
        loc1 = get_location(st1)
        loc2 = get_location(st2)
        if loc1 == loc2:
            return True
        elif isinstance(loc1, list):
            if set(loc2).issubset(set(loc1)):
                return True
    return False


def location_refinement(st1, st2, hierarchies):
    """Return True if there is a location-aware refinement between stmts."""
    if type(st1) != type(st2):
        return False
    if isinstance(st1, Event):
        event_ref = event_location_refinement(st1, st2, hierarchies)
        return event_ref
    elif isinstance(st1, Influence):
        subj_ref = event_location_refinement(st1.subj, st2.subj,
                                             hierarchies)
        obj_ref = event_location_refinement(st1.obj, st2.obj,
                                            hierarchies)
        return subj_ref and obj_ref
    else:
        return st1.refinement_of(st2, hierarchies)


def event_location_time_matches(event):
    """Return Event matches key which takes location and time into account."""
    mk = location_matches(event)
    if not has_time(event):
        return mk
    time = get_time(event)
    matches_key = str((mk, time.start, time.end, time.duration))
    return matches_key


def location_time_matches(stmt):
    """Return matches key which takes location and time into account."""
    if isinstance(stmt, Event):
        return event_location_time_matches(stmt)
    elif isinstance(stmt, Influence):
        subj_mk = event_location_time_matches(stmt.subj)
        obj_mk = event_location_time_matches(stmt.obj)
        return str((stmt.matches_key(), subj_mk, obj_mk))
    else:
        return stmt.matches_key()


def event_location_time_refinement(st1, st2, hierarchies):
    """Return True if there is a location/time refinement between Events."""
    ref = location_refinement(st1, st2, hierarchies)
    if not ref:
        return False
    if not has_time(st2):
        return True
    elif not has_time(st1):
        return False
    else:
        return st1.context.time.refinement_of(st2.context.time)


def location_time_refinement(st1, st2, hierarchies):
    """Return True if there is a location/time refinement between stmts."""
    if type(st1) != type(st2):
        return False
    if isinstance(st1, Event):
        return event_location_time_refinement(st1, st2, hierarchies)
    elif isinstance(st1, Influence):
        ref = st1.refinement_of(st2, hierarchies)
        if not ref:
            return False
        subj_ref = event_location_time_refinement(st1.subj, st2.subj,
                                                  hierarchies)
        obj_ref = event_location_time_refinement(st1.obj, st2.obj,
                                                 hierarchies)
        return subj_ref and obj_ref


def agent_grounding_matches(agent):
    """Return an Agent matches key just based on grounding, not state."""
    if agent is None:
        return None
    return str(agent.entity_matches_key())


def agents_stmt_type_matches(stmt):
    """Return a matches key just based on Agent grounding and Stmt type."""
    agents = [agent_grounding_matches(a) for a in stmt.agent_list()]
    key = str((stmt.__class__.__name__, agents))
    return key


def has_delta(stmt):
    if not stmt.delta:
        return False
    return True


def get_delta(stmt):
    delta = stmt.delta
    if isinstance(delta, QualitativeDelta):
        return delta.polarity
    elif isinstance(delta, QuantitativeState):
        return (delta.entity, delta.value, delta.unit, delta.polarity)


def event_location_time_delta_matches(event):
    mk = event_location_time_matches(event)
    if not has_delta:
        return mk
    delta = get_delta(event)
    matches_key = str((mk, delta))
    return matches_key


def location_time_delta_matches(stmt):
    if isinstance(stmt, Event):
        return event_location_time_delta_matches(stmt)
    elif isinstance(stmt, Influence):
        subj_mk = event_location_time_delta_matches(stmt.subj)
        obj_mk = event_location_time_delta_matches(stmt.obj)
        return str((stmt.matches_key(), subj_mk, obj_mk))
    else:
        return stmt.matches_key()


def event_location_time_delta_refinement(st1, st2, hierarchies):
    loc_time_ref = event_location_time_refinement(st1, st2, hierarchies)
    if not loc_time_ref:
        return False
    if not st2.delta:
        return True
    elif not st1.delta:
        return False
    else:
        return st1.delta.refinement_of(st2.delta)


def location_time_delta_refinement(st1, st2, hierarchies):
    if isinstance(st1, Event):
        return event_location_time_delta_refinement(st1, st2, hierarchies)
    elif isinstance(st1, Influence):
        ref = st1.refinement_of(st2, hierarchies)
        if not ref:
            return False
        subj_ref = event_location_time_delta_refinement(st1.subj, st2.subj,
                                                        hierarchies)
        obj_ref = event_location_time_delta_refinement(st1.obj, st2.obj,
                                                       hierarchies)
        return subj_ref and obj_ref
    else:
        return st1.refinement_of(st2, hierarchies)
