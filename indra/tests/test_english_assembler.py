import indra.assemblers.english.assembler as ea
from indra.statements import *


def _sustring_by_coords(text, coords):
    return text[coords[0]:coords[1]]


def test_agent_basic():
    ag = ea._assemble_agent_str(Agent('EGFR'))
    assert isinstance(ag, ea.AgentWithCoordinates)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mod():
    mc = ModCondition('phosphorylation')
    a = Agent('EGFR', mods=mc)
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'phosphorylated EGFR'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mod2():
    mc = ModCondition('phosphorylation', 'tyrosine')
    a = Agent('EGFR', mods=mc)
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'tyrosine-phosphorylated EGFR'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mod3():
    mc = ModCondition('phosphorylation', 'tyrosine', '1111')
    a = Agent('EGFR', mods=mc)
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR phosphorylated on Y1111'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mods():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation', 'tyrosine', '1234')
    a = Agent('EGFR', mods=[mc1, mc2])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR phosphorylated on Y1111 and Y1234'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mods2():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation', 'tyrosine')
    a = Agent('EGFR', mods=[mc1, mc2])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR phosphorylated on Y1111 and tyrosine'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_mods3():
    mc1 = ModCondition('phosphorylation', 'tyrosine', '1111')
    mc2 = ModCondition('phosphorylation')
    a = Agent('EGFR', mods=[mc1, mc2])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR phosphorylated on Y1111 and an unknown residue'
    assert ag.agent_str[ag.coords[0]:ag.coords[1]] == 'EGFR'


def test_agent_mods_pos_only():
    mc1 = ModCondition('phosphorylation', None, '1111')
    a = Agent('EGFR', mods=[mc1])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR phosphorylated on amino acid 1111'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_bound():
    bc = BoundCondition(Agent('EGF'), True)
    a = Agent('EGFR', bound_conditions=[bc])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR bound to EGF'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_not_bound():
    bc = BoundCondition(Agent('EGF'), False)
    a = Agent('EGFR', bound_conditions=[bc])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR not bound to EGF'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_bound_two():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), True)
    a = Agent('EGFR', bound_conditions=[bc, bc2])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR bound to EGF and EGFR'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_bound_three():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), True)
    bc3 = BoundCondition(Agent('GRB2'), True)
    a = Agent('EGFR', bound_conditions=[bc, bc2, bc3])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR bound to EGF, EGFR, and GRB2'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_agent_bound_mixed():
    bc = BoundCondition(Agent('EGF'), True)
    bc2 = BoundCondition(Agent('EGFR'), False)
    a = Agent('EGFR', bound_conditions=[bc, bc2])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'EGFR bound to EGF and not bound to EGFR'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'EGFR'


def test_phos_noenz():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a)
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'MAP2K1 is phosphorylated.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'MAP2K1'


def test_phos_noenz2():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a, 'serine')
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'MAP2K1 is phosphorylated on serine.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'MAP2K1'


def test_phos_noenz3():
    a = Agent('MAP2K1')
    st = Phosphorylation(None, a, 'serine', '222')
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'MAP2K1 is phosphorylated on S222.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'MAP2K1'


def test_phos_enz():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    st = Phosphorylation(b, a, 'serine', '222')
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'BRAF phosphorylates MAP2K1 on S222.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'BRAF'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'MAP2K1'


def test_phos_indirect():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    ev = Evidence(epistemics={'direct': False})
    st = Phosphorylation(b, a, 'serine', '222', evidence=[ev])
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'BRAF leads to the phosphorylation of MAP2K1 on S222.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'BRAF'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'MAP2K1'


def test_phos_enz2():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st = Dephosphorylation(b, a, 'serine', '222')
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'PP2A dephosphorylates MAP2K1 on S222.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'PP2A'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'MAP2K1'


def test_dephos_no_residue():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st = Dephosphorylation(b, a, None, '222')
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'PP2A dephosphorylates MAP2K1 at position 222.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'PP2A'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'MAP2K1'


def test_ubiq_stmt():
    st = Ubiquitination(Agent('X'), Agent('Y'))
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'X ubiquitinates Y.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'X'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'Y'


def test_deubiq_stmt():
    st = Deubiquitination(Agent('X'), Agent('Y'))
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'X deubiquitinates Y.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'X'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'Y'


def test_deubiq_noenz():
    st = Deubiquitination(None, Agent('Y'))
    sb = ea._assemble_modification(st)
    print(sb.sentence)
    assert sb.sentence == 'Y is deubiquitinated.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'Y'


def test_complex_one():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    st = Complex([a, b])
    sb = ea._assemble_complex(st)
    print(sb.sentence)
    assert sb.sentence == 'MAP2K1 binds BRAF.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'MAP2K1'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'BRAF'


def test_complex_more():
    a = Agent('MAP2K1')
    b = Agent('BRAF')
    c = Agent('RAF1')
    st = Complex([a, b, c])
    sb = ea._assemble_complex(st)
    print(sb.sentence)
    assert sb.sentence == 'MAP2K1 binds BRAF and RAF1.'
    assert _sustring_by_coords(sb.sentence, sb.agents[0].coords) == 'MAP2K1'
    assert _sustring_by_coords(sb.sentence, sb.agents[1].coords) == 'BRAF'
    assert _sustring_by_coords(sb.sentence, sb.agents[2].coords) == 'RAF1'


def test_assemble_one():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st = Dephosphorylation(b, a, 'serine', 222)
    s = _stmt_to_text(st)
    assert s == 'PP2A dephosphorylates MAP2K1 on S222.'


def test_assemble_more():
    a = Agent('MAP2K1')
    b = Agent('PP2A')
    st1 = Dephosphorylation(b, a, 'serine', 222)
    b = Agent('BRAF')
    c = Agent('RAF1')
    st2 = Complex([a, b, c])
    e = ea.EnglishAssembler()
    e.add_statements([st1, st2])
    s = e.make_model()
    print(s)
    assert s == \
        'PP2A dephosphorylates MAP2K1 on S222. MAP2K1 binds BRAF and RAF1.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'PP2A'
    assert _sustring_by_coords(s, e.stmt_agents[0][1].coords) == 'MAP2K1'
    assert _sustring_by_coords(s, e.stmt_agents[1][0].coords) == 'MAP2K1'
    assert _sustring_by_coords(s, e.stmt_agents[1][1].coords) == 'BRAF'
    assert _sustring_by_coords(s, e.stmt_agents[1][2].coords) == 'RAF1'


def test_autophos():
    a = Agent('EGFR')
    st = Autophosphorylation(a, 'Y')
    s = _stmt_to_text(st)
    assert s == 'EGFR phosphorylates itself on tyrosine.'


def test_activation():
    st = Activation(Agent('MEK'), Agent('ERK'))
    s = _stmt_to_text(st)
    assert s == 'MEK activates ERK.'


def test_inhibition():
    st = Inhibition(Agent('MEK'), Agent('ERK'))
    s = _stmt_to_text(st)
    assert s == 'MEK inhibits ERK.'


def test_regulateamount():
    st = IncreaseAmount(Agent('TP53'), Agent('MDM2'))
    s = _stmt_to_text(st)
    assert s == 'TP53 increases the amount of MDM2.'
    st = DecreaseAmount(Agent('TP53'), Agent('MDM2'))
    s = _stmt_to_text(st)
    assert s == 'TP53 decreases the amount of MDM2.'
    st = DecreaseAmount(None, Agent('MDM2'))
    s = _stmt_to_text(st)
    assert s == 'MDM2 is degraded.'
    st = IncreaseAmount(None, Agent('MDM2'))
    s = _stmt_to_text(st)
    assert s == 'MDM2 is produced.'


def test_agent_loc():
    a = Agent('BRAF', location='cytoplasm')
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'BRAF in the cytoplasm'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'BRAF'


def test_agent_mut():
    a = Agent('BRAF', mutations=[MutCondition('600','V', 'E')])
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'BRAF-V600E'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'BRAF'


def test_agent_mut_plus():
    a = Agent('BRAF', mutations=[MutCondition('600','V', 'E')],
              location='nucleus')
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'BRAF-V600E in the nucleus'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'BRAF'


def test_agent_activity():
    a = Agent('BRAF', activity=ActivityCondition('activity', True))
    ag = ea._assemble_agent_str(a)
    print(ag.agent_str)
    assert ag.agent_str == 'active BRAF'
    assert _sustring_by_coords(ag.agent_str, ag.coords) == 'BRAF'


def test_agent_activity_stmt():
    braf = Agent('BRAF', activity=ActivityCondition('activity', True))
    mek = Agent('MAP2K1')
    st = Activation(braf, mek)
    e = ea.EnglishAssembler()
    e.add_statements([st])
    s = e.make_model()
    print(s)
    assert s == 'Active BRAF activates MAP2K1.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'BRAF'
    assert _sustring_by_coords(s, e.stmt_agents[0][1].coords) == 'MAP2K1'


def test_translocation():
    st1 = Translocation(Agent('FOXO3A'))
    st2 = Translocation(Agent('FOXO3A'), 'cytoplasm')
    st3 = Translocation(Agent('FOXO3A'), None, 'nucleus')
    st4 = Translocation(Agent('FOXO3A'), 'cytoplasm', 'nucleus')
    e = ea.EnglishAssembler()
    e.add_statements([st1])
    s = e.make_model()
    assert s == 'FOXO3A translocates.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'FOXO3A'
    e = ea.EnglishAssembler()
    e.add_statements([st2])
    s = e.make_model()
    assert s == 'FOXO3A translocates from the cytoplasm.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'FOXO3A'
    e = ea.EnglishAssembler()
    e.add_statements([st3])
    s = e.make_model()
    assert s == 'FOXO3A translocates to the nucleus.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'FOXO3A'
    e = ea.EnglishAssembler()
    e.add_statements([st4])
    s = e.make_model()
    assert s == 'FOXO3A translocates from the cytoplasm to the nucleus.'
    assert _sustring_by_coords(s, e.stmt_agents[0][0].coords) == 'FOXO3A'


def test_gef():
    st = Gef(Agent('SOS1'), Agent('KRAS'))
    s = _stmt_to_text(st)
    assert s == 'SOS1 is a GEF for KRAS.'


def test_gap():
    st = Gap(Agent('RASA1'), Agent('KRAS'))
    s = _stmt_to_text(st)
    assert s == 'RASA1 is a GAP for KRAS.'


def test_methylation():
    st = Methylation(None, Agent('SLF11'))
    s = _stmt_to_text(st)
    assert s == 'SLF11 is methylated.'


def test_generic_mod_state():
    mc = ModCondition('modification')
    st = Activation(Agent('MEK', mods=[mc]), Agent('ERK'))
    s = _stmt_to_text(st)
    assert s == 'Modified MEK activates ERK.'


def test_generic_mutation():
    mc = MutCondition(None, None, None)
    st = Activation(Agent('MEK', mutations=[mc]), Agent('ERK'))
    s = _stmt_to_text(st)
    assert s == 'Mutated MEK activates ERK.'


def test_activity_conditions():
    def to_english(act_type, is_act):
        ac = ActivityCondition(act_type, is_act)
        st = Activation(Agent('MEK', activity=ac), Agent('ERK'))
        s = _stmt_to_text(st)
        return s

    assert to_english('activity', True) == 'Active MEK activates ERK.'
    assert to_english('activity', False) == 'Inactive MEK activates ERK.'
    assert to_english('gtpbound', True) == \
        'GTP-bound active MEK activates ERK.'
    assert to_english('gtpbound', False) == \
        'GDP-bound inactive MEK activates ERK.'
    assert to_english('catalytic', False) == \
        'Catalytically inactive MEK activates ERK.'
    assert to_english('kinase', True) == \
        'Kinase-active MEK activates ERK.'
    assert to_english('gap', True) == \
        'GAP-active MEK activates ERK.'
    assert to_english('gef', False) == \
        'GEF-inactive MEK activates ERK.'


def test_conversion():
    st = Conversion(Agent('RAS'), [Agent('GTP')], [Agent('GDP')])
    s = _stmt_to_text(st)
    assert s == 'RAS catalyzes the conversion of GTP into GDP.'

    st = Conversion(Agent('RAS'), [Agent('GTP'), Agent('X')],
                    [Agent('GDP'), Agent('Y')])
    s = _stmt_to_text(st)
    assert s == 'RAS catalyzes the conversion of GTP and X into GDP and Y.'

    st = Conversion(None, [Agent('GTP'), Agent('X')], [Agent('GDP'), Agent('Y')])
    s = _stmt_to_text(st)
    assert s == 'GTP and X is converted into GDP and Y.'


def test_influence():
    st = Influence(Event(Concept('food')), Event(Concept('hunger'),
                                                 delta=QualitativeDelta(
                                                     polarity=-1)))
    s = _stmt_to_text(st)
    assert s == 'Food causes a decrease in hunger.', s

    st = Influence(Event(Concept('food'), delta=QualitativeDelta(
                                            polarity=1)),
                   Event(Concept('hunger'), delta=QualitativeDelta(
                                            polarity=-1)))
    s = _stmt_to_text(st)
    assert s == 'An increase in food causes a decrease in hunger.', s

    st = Influence(Event(Concept('food'), delta=QualitativeDelta(
                                            polarity=-1)),
                   Event(Concept('hunger'), delta=QualitativeDelta(
                                            polarity=1)))
    s = _stmt_to_text(st)
    assert s == 'A decrease in food causes an increase in hunger.', s


def test_association():
    st = Association([Event(Concept('food')), Event(Concept('hunger'))])
    s = _stmt_to_text(st)
    assert s == 'Food is associated with hunger.'


def test_active_form():
    st = ActiveForm(Agent('BRAF',
                          mods=[ModCondition('phosphorylation', 'T', '396')]),
                    'kinase', True)
    s = _stmt_to_text(st)
    assert s == 'BRAF phosphorylated on T396 is kinase-active.'


def test_get_base_verb():
    assert ea.statement_base_verb('inhibition') == 'inhibit'
    assert ea.statement_base_verb('dephosphorylation') == 'dephosphorylate'
    assert ea.statement_base_verb('complex') == 'bind'


def _stmt_to_text(st):
    e = ea.EnglishAssembler()
    e.add_statements([st])
    s = e.make_model()
    print(s)
    return s
