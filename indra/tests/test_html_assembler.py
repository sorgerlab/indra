import re
from indra.statements import *
from indra.assemblers.english import AgentWithCoordinates
from indra.assemblers.html.assembler import HtmlAssembler, tag_text, loader, \
    _format_evidence_text, tag_agents
from indra.util.statement_presentation import AveAggregator, StmtStat, StmtGroup


def make_stmt():
    src = Agent('SRC', db_refs={'HGNC': '11283'})
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    ev = Evidence(text="We noticed that the Src kinase was able to "
                       "phosphorylate Ras proteins.",
                  source_api='test', pmid='1234567',
                  annotations={'agents': {'raw_text': ['Src kinase',
                                                       'Ras proteins']}})
    st = Phosphorylation(src, ras, 'tyrosine', '32', evidence=[ev])
    return st


def make_bad_stmt():
    subj = None  # None agent
    # Agent without name and several matches in db_refs
    ras = Agent('', db_refs={'FPLX': {'RAS', 'Ras'}, 'TEXT': 'RAS'})
    ev = Evidence(text="Ras is phosphorylated",
                  source_api='test', pmid='1234',
                  annotations={'agents': {'raw_text': [None, None]}})  # no raw
    st = Phosphorylation(subj, ras, 'tyrosine', '32', evidence=[ev])
    return st


def test_format_evidence_text():
    stmt = make_stmt()
    ev_list = _format_evidence_text(stmt)
    assert len(ev_list) == 1
    ev = ev_list[0]
    assert isinstance(ev, dict)
    assert set(ev.keys()) == {'source_api', 'text_refs', 'text', 'source_hash',
                              'pmid', 'num_curations', 'num_correct',
                              'num_incorrect'}
    assert ev['source_api'] == 'test'
    assert ev['text_refs']['PMID'] == '1234567'
    assert ev['text'] == ('We noticed that the '
                          '<span class="badge badge-subject">Src kinase</span> '
                          'was able to phosphorylate '
                          '<span class="badge badge-object">'
                          'Ras proteins</span>.'), ev['text']


def test_assembler():
    stmt = make_stmt()
    ha = HtmlAssembler([stmt])
    result = ha.make_model()
    assert isinstance(result, str)
    # Read from the template file and make sure the beginning and end of the
    # content matches
    template, _, _ = loader.get_source(None, 'indra/template.html')
    assert result.startswith(template[0:100])
    # Make sure assembler works with other parameters provided
    stmt2 = make_bad_stmt()
    ha = HtmlAssembler(
        source_counts={stmt.get_hash(): {'test': 1},
                       stmt2.get_hash(): {'test': 1}},
        ev_counts={stmt.get_hash(): 1, stmt2.get_hash(): 1},
        db_rest_url='test.db.url')
    ha.add_statements([stmt, stmt2])
    result = ha.make_model(grouping_level='agent-pair')
    assert isinstance(result, str)
    result = ha.make_model(grouping_level='statement')
    assert isinstance(result, str)
    # Make sure warning can be appended
    ha.append_warning('warning')
    assert ('\t<span style="color:red;">(CAUTION: warning occurred when '
            'creating this page.)</span>' in ha.model)
    # Make sure model is created before saving
    ha = HtmlAssembler([stmt])
    assert not ha.model
    ha.save_model('tempfile.html')
    assert ha.model


def test_tag_agents():
    # tag_agents input can be either regular agents or agents with coordinates
    english = 'SRC phosphorylates RAS on Y32.'
    src = Agent('SRC', db_refs={'HGNC': '11283'})
    ras = Agent('RAS', db_refs={'FPLX': 'RAS'})
    src_with_coords = AgentWithCoordinates(
        'SRC', 'SRC', db_refs={'HGNC': '11283'}, coords=(0, 3))
    ras_with_coords = AgentWithCoordinates(
        'RAS', 'RAS', db_refs={'FPLX': 'RAS'}, coords=(19, 22))
    assert tag_agents(english, [src, ras]) == tag_agents(
        english, [src_with_coords, ras_with_coords])


def test_tag_text():
    """If there are overlapping or nested matches, show only one."""
    text = 'FooBarBaz binds Foo.'
    indices = []
    for span in ('FooBarBaz', 'Foo'):
        tag_start = "<%s>" % span
        tag_close = "</%s>" % span
        indices += [(m.start(), m.start() + len(span), span,
                     tag_start, tag_close)
                    for m in re.finditer(re.escape(span), text)]
    tagged_text = tag_text(text, indices)
    assert tagged_text == '<FooBarBaz>FooBarBaz</FooBarBaz> binds ' \
                          '<Foo>Foo</Foo>.'


def test_tag_bad_text():
    ev = Evidence("bogus", text="<Foo> binds Bar& (<10 & >20)",
                  annotations={"agents": {"raw_text": ["<Foo>", "Bar&"]}})
    stmt = Complex([Agent("Foo"), Agent("Bar")], evidence=[ev])
    ev_list = _format_evidence_text(stmt)
    fmt_ev = ev_list[0]
    assert fmt_ev['text'] == ("<span class=\"badge badge-other\">&lt;Foo&gt;"
                              "</span> binds <span class=\"badge badge-other\">"
                              "Bar&amp;</span> (&lt;10 &amp; &gt;20)")


def test_influence():
    c2 = Concept('food insecurity',
                 db_refs={'WM': [('wm/food_insecurity', 1.0)]})
    c1 = Concept('floods',
                 db_refs={'WM': [('wm/floods', 1.0)]})
    c3 = Concept('x', db_refs={})
    stmt = Influence(Event(c1), Event(c2))
    stmt2 = Influence(Event(c1), Event(c3))
    ha = HtmlAssembler([stmt, stmt2])
    ha.make_model()


def test_active_form():
    stmt = ActiveForm(Agent('MAPK1', mods=[ModCondition('phosphorylation')]),
                      'kinase', True)
    ha = HtmlAssembler([stmt])
    ha.make_model()
    # Case when it's not active
    stmt = ActiveForm(Agent('MAPK1', mods=[ModCondition('phosphorylation')]),
                      'activity', False)
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_complex():
    stmt = Complex([Agent('BRAF'), Agent('RAF1')])
    ha = HtmlAssembler([stmt])
    ha.make_model()
    # Complex with more than two members
    stmt = Complex([Agent('BRAF'), Agent('RAF1'), Agent('YWAH')])
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_conversion():
    stmt = Conversion(Agent('RAS'), [Agent('GTP'), Agent('B')],
                      [Agent('GDP'), Agent('C')])
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_has_activity():
    stmt = HasActivity(Agent('MAPK1'), 'activity', True)
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_phosphorylation():
    stmt = Phosphorylation(Agent('MAP2K1'), Agent('MAPK1'), 'T', '185')
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_autophosphorylation():
    stmt = Autophosphorylation(
        Agent('P38', bound_conditions=[BoundCondition(Agent('TAB1'))]))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_dephosphorylation():
    stmt = Dephosphorylation(Agent('DUSP6'), Agent('MAPK1'), 'T', '185')
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_inhibition():
    stmt = Inhibition(Agent('DUSP4'), Agent('MAPK1'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_activation():
    stmt = Activation(Agent('MAP2K1'), Agent('MAPK1'), 'kinase')
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_gef():
    stmt = Gef(Agent('SOS1'), Agent('KRAS'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_gap():
    stmt = Gap(Agent('RASA1'), Agent('KRAS'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_translocation():
    stmt = Translocation(Agent('FOXO3A'), None, 'nucleus')
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_increase_amount():
    stmt = IncreaseAmount(Agent('TP53'), Agent('MDM2'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_decrease_amount():
    stmt = DecreaseAmount(Agent('TP53'), Agent('MDM2'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_association():
    stmt = Association([Event(Concept('a')), Event(Concept('b'))])
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_event():
    stmt = Event(Concept('a'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def test_migration():
    stmt = Migration(Concept('migration'))
    ha = HtmlAssembler([stmt])
    ha.make_model()


def _get_sort_corpus():
    # Create a list of statements with bogus agents.
    stmts = [
        Inhibition(
            Agent('Fez'), Agent('Baz'),
            evidence=[Evidence('medscan', text="Fez-|Baz"),
                      Evidence('medscan', text="Baz|-Fez"),
                      Evidence('reach', text="Fez inhibits Baz."),
                      Evidence('isi', text="Baz is inhibited by Fez.")]
        ),
        DecreaseAmount(
            Agent('Fez'), Agent('Baz'),
            evidence=[Evidence('reach', text="Fez->Baz")]
        ),
        Complex(
            [Agent('Fez'), Agent('Baz'), Agent('Bar')],
            evidence=[Evidence('sparser', text="Fez-Baz-Bar complex"),
                      Evidence('reach', text="Complex of Fez, Baz, & Bar")]
        ),
        Phosphorylation(
            Agent('Bar'), Agent('Baz'), 'T', '185',
            evidence=[Evidence('reach',
                               text="Bar phosphorylates Baz on T 46"),
                      Evidence('trips',
                               text="Bar phosphorylate Baz on Tyrosine "
                                    "forty-six")]
        ),
        Phosphorylation(
            Agent('Bar'), Agent('Baz'),
            evidence=[Evidence('sparser',
                               text="Bar phosphorylates Baz"),
                      Evidence('isi',
                               text="Bar phosphorylates Baz."),
                      Evidence('sparser',
                               text="Baz is phosphorylated by Bar.")]
        ),
        Conversion(Agent('Fez'), [Agent('Far'), Agent('Faz')],
                   [Agent('Bar'), Agent('Baz')],
                   evidence=[Evidence('reach',
                                      text='Fez converts Far and Faz into Bar '
                                           'and Baz.')])
    ]

    # Set belief values.
    beliefs = [0.9, 0.8, 0.6, 0.99, 1, 0.75]
    for stmt, belief in zip(stmts, beliefs):
        stmt.belief = belief
    return stmts


def test_sort_default():
    ha = HtmlAssembler(_get_sort_corpus())

    # Test the ordering of the statements in the default mode of make_json_model
    json_model = ha.make_json_model()
    assert list(json_model.keys()) == ['Fez-Baz', 'Bar-Baz', 'Fez-Bar']
    exp_stmt_counts = {'Fez-Baz': 4, 'Bar-Baz': 2, 'Fez-Bar': 2}
    assert all(len(json_model[k]['stmts_formatted']) == n
               for k, n in exp_stmt_counts.items())
    ev_counts = {k: sum(len(s['evidence']) for r in m['stmts_formatted']
                        for s in r['stmt_info_list'])
                 for k, m in json_model.items()}
    assert ev_counts == {'Fez-Baz': 8, 'Bar-Baz': 7, 'Fez-Bar': 3}, ev_counts

    # Check to make sure the HTML assembler runs.
    model = ha.make_model()
    with open('test_agent_pair.html', 'w') as f:
        f.write(model)


def test_sort_group_by_relation():
    ha = HtmlAssembler(_get_sort_corpus())

    # Test ordering, grouping by relation.
    json_model = ha.make_json_model(grouping_level='relation')
    assert list(json_model.keys()) == ['all-relations']
    relations = json_model['all-relations']['stmts_formatted']
    assert len(relations) == 5, len(relations)

    # Make sure the HTML assembles.
    model = ha.make_model(grouping_level='relation')
    with open('test_relation.html', 'w') as f:
        f.write(model)


def test_sort_group_by_statement():
    ha = HtmlAssembler(_get_sort_corpus())

    # Test ordering and grouping by statement.
    json_model = ha.make_json_model(grouping_level='statement')
    assert list(json_model.keys()) == ['all-statements']
    assert len(json_model['all-statements']['stmts_formatted']) == 1
    statements = \
        json_model['all-statements']['stmts_formatted'][0]['stmt_info_list']
    assert len(statements) == 6
    assert [len(s['evidence']) for s in statements] == [4, 3, 2, 2, 1, 1]

    # Make sure the html assembly works.
    ha.make_model(grouping_level='statement')


def test_sort_group_by_statement_sort_by_none():
    stmts = _get_sort_corpus()
    ha = HtmlAssembler(stmts, sort_by=None)

    json_model = ha.make_json_model(grouping_level='statement')
    statements = \
        json_model['all-statements']['stmts_formatted'][0]['stmt_info_list']
    got_h_list = [int(s['hash']) for s in statements]
    inp_h_list = [s.get_hash() for s in stmts]
    assert got_h_list == inp_h_list


def test_sort_group_by_statement_custom_ordering():
    stmts = _get_sort_corpus()

    custom_values = [0.1, 0.2, 0.15, 0.6, 0.3, 0.8]
    val_dict = {s.get_hash(): v for v, s in zip(custom_values, stmts)}

    custom_stat = StmtStat('value', val_dict, float, AveAggregator)
    stmt_gather = StmtGroup.from_stmt_stats(custom_stat)

    ha = HtmlAssembler(stmts, sort_by='value', stmt_stat_gather=stmt_gather)
    json_model = ha.make_json_model(grouping_level='statement')

    statements = \
        json_model['all-statements']['stmts_formatted'][0]['stmt_info_list']
    got_h_list = [int(s['hash']) for s in statements]
    exp_h_list = sorted((h for h in val_dict.keys()), key=lambda h: val_dict[h],
                        reverse=True)
    assert got_h_list == exp_h_list

    ha.make_model(grouping_level='statement')


def test_sort_group_by_relation_custom_ordering():
    stmts = _get_sort_corpus()

    custom_values = [0.1, 0.2, 0.15, 0.6, 0.3, 0.8]
    val_dict = {s.get_hash(): v for v, s in zip(custom_values, stmts)}

    custom_stat = StmtStat('value', val_dict, float, AveAggregator)
    stmt_gather = StmtGroup.from_stmt_stats(custom_stat)

    ha = HtmlAssembler(stmts, sort_by='value', stmt_stat_gather=stmt_gather)
    json_model = ha.make_json_model(grouping_level='relation')
    assert list(json_model.keys()) == ['all-relations']
    relations = json_model['all-relations']['stmts_formatted']
    assert len(relations) == 5, len(relations)
    relation_names = [rel['short_name'] for rel in relations]
    exp_relation_names = [
        '<b>Fez</b> catalyzes the conversion of <b>Far</b> and <b>Faz</b> into '
        '<b>Bar</b> and <b>Baz</b>.',
        '<b>Bar</b> phosphorylates <b>Baz</b>.',
        '<b>Fez</b> decreases the amount of <b>Baz</b>.',
        '<b>Bar</b> binds <b>Baz</b> and <b>Fez</b>.',
        '<b>Fez</b> inhibits <b>Baz</b>.'
    ]
    assert relation_names == exp_relation_names

    ha.make_model(grouping_level='relation')


def test_sort_group_by_agent_custom_ordering():
    stmts = _get_sort_corpus()

    custom_values = [0.1, 0.2, 0.15, 0.6, 0.3, 0.8]
    val_dict = {s.get_hash(): v for v, s in zip(custom_values, stmts)}

    custom_stat = StmtStat('value', val_dict, float, AveAggregator)
    stmt_gather = StmtGroup.from_stmt_stats(custom_stat)

    ha = HtmlAssembler(stmts, sort_by='value', stmt_stat_gather=stmt_gather)
    json_model = ha.make_json_model(grouping_level='agent-pair')
    assert len(json_model.keys()) == 4

    # This result was slightly counter-intuitive, but recall that averages will
    # mean a grouping with the conversion will always have a lower value than
    # the conversion itself, so it makes sense for it to come out on top.
    assert list(json_model.keys()) == ['Fez-Far-Faz-Bar-Baz', 'Fez-Bar',
                                       'Bar-Baz', 'Fez-Baz']

    ha.make_model(grouping_level='agent-pair')


def test_sort_group_by_statement_custom_function():
    stmts = _get_sort_corpus()

    ha = HtmlAssembler(stmts,
                       sort_by=lambda d: 4*d['trips'] + 2*d['reach']
                                         + 2*d['medscan'] + d['sparser']
                                         - d['isi'])
    json_model = ha.make_json_model(grouping_level='statement')
    statements = \
        json_model['all-statements']['stmts_formatted'][0]['stmt_info_list']
    assert len(statements) == len(stmts)
    exp_order = ['6106301533612997', '-17995265549545446', '34182032179844940',
                 '32266861591785935', '-30059881887512900', '-5998595995539618']
    assert [s['hash'] for s in statements] == exp_order

    ha.make_model(grouping_level='statement')


def test_sort_group_by_relation_custom_function():
    stmts = _get_sort_corpus()

    ha = HtmlAssembler(stmts,
                       sort_by=lambda d: 4*d['trips'] + 2*d['reach']
                                         + 2*d['medscan'] + d['sparser']
                                         - d['isi'])
    json_model = ha.make_json_model(grouping_level='relation')
    relations = json_model['all-relations']['stmts_formatted']
    assert len(relations) == 5, len(relations)
    relation_names = [rel['short_name'] for rel in relations]
    exp_rel_names = [
        '<b>Bar</b> phosphorylates <b>Baz</b>.',
        '<b>Fez</b> inhibits <b>Baz</b>.',
        '<b>Bar</b> binds <b>Baz</b> and <b>Fez</b>.',
        '<b>Fez</b> decreases the amount of <b>Baz</b>.',
        '<b>Fez</b> catalyzes the conversion of <b>Far</b> and <b>Faz</b> into '
        '<b>Bar</b> and <b>Baz</b>.'
    ]
    assert relation_names == exp_rel_names, relation_names

    ha.make_model(grouping_level='relation')


def test_sort_group_by_agent_pair_custom_function():
    stmts = _get_sort_corpus()

    ha = HtmlAssembler(stmts,
                       sort_by=lambda d: 4*d['trips'] + 2*d['reach']
                                         + 2*d['medscan'] + d['sparser']
                                         - d['isi'])
    json_model = ha.make_json_model(grouping_level='agent-pair')
    assert list(json_model.keys()) == ['Fez-Baz', 'Bar-Baz', 'Fez-Bar']

    ha.make_model(grouping_level='agent-pair')
