
def find_rules(model, monomer_name):
    rules = []
    for r in model.rules:
        lhs_cps = r.reactant_pattern.complex_patterns
        rhs_cps = r.product_pattern.complex_patterns
        found = False
        for cp in lhs_cps + rhs_cps:
            for mp in cp.monomer_patterns:
                if mp.monomer.name == monomer_name:
                    found = True
        if found:
            rules.append(r)

    return rules
