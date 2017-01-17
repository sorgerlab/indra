import pickle
import logging
from copy import deepcopy
from indra.statements import *

logger = logging.getLogger('active_forms')

def apply_active_forms(stmts):
    # Collect Modification statements where the enzyme has no conditions
    active_forms = {}
    blank_stmts = []
    other_stmts = []
    for stmt in stmts:
        if isinstance(stmt, Modification):
            enz = stmt.enz
            if enz is not None and \
                not (enz.mods or enz.bound_conditions or enz.mutations or
                                                 enz.location or enz.active):
                blank_stmts.append(stmt)
                logger.info("Blank: %s" % stmt)
            else:
                other_stmts.append(stmt)
        elif isinstance(stmt, ActiveForm):
            active_form_list = active_forms.get(stmt.agent.name, [])
            active_form_list.append(stmt)
            active_forms[stmt.agent.name] = active_form_list
        else:
            other_stmts.append(stmt)

    new_stmts = []
    for stmt in blank_stmts:
        enz_name = stmt.enz.name
        active_form_list = active_forms.get(enz_name)
        # No active form info for this enzyme, pass it through unchanged
        if not active_form_list:
            new_stmts.append(stmt)
        # We've got some different active forms, make new stmts for each one
        else:
            # Copy the agent in the active form
            for af in active_form_list:
                agent_copy = deepcopy(af.agent)
                mod_copy = deepcopy(stmt)
                mod_copy.enz = agent_copy
                new_stmts.append(mod_copy)

    # Filters out all active forms?
    all_stmts = other_stmts + new_stmts + list(active_forms.values())
    return all_stmts

if __name__ == '__main__':
    with open(sys.argv[1], 'rb') as f:
        stmts = pickle.load(f)
    for s in stmts:
        logger.info("Before: %s" % s)
    new_stmts = apply_active_forms(stmts)
    for s in new_stmts:
        logger.info("After: %s" % s)

    print("Before %d, after %d" % (len(stmts), len(new_stmts)))
