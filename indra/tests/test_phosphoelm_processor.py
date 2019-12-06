from indra.statements import Phosphorylation
from indra.sources.phosphoelm.processor import PhosphoElmProcessor, \
    _agent_from_str
from indra.sources.phosphoelm.api import _get_json_from_entry_rows

columns = ['acc', 'sequence', 'position', 'code', 'pmids', 'kinases',
           'source', 'species', 'entry_date']
non_human_no_kinase = ['O08539',
                       'FAKEGSKG', '6', 'S', '17114649', '', 'HTP',
                       'Mus musculus', '2005-03-14 12:16:11.108314+01']
human_no_kinase = ['O14543',  # SOCS3
                   'FAKESEQUENCESRPLDTSLRLKTFSSKSEYQL', '31', 'Y',
                   '12783885', '', 'LTP', 'Homo sapiens',
                   '2006-10-17 12:06:48.271076+02']
human_kinase1 = ['O14543',  # SOCS3
                 'FAKESEQUENCESRPLDTSLRLKTFSSKSEYQL',
                 '31', 'Y', '12783885', 'Lck', 'LTP', 'Homo sapiens',
                 '2006-10-17 12:06:48.16767+02']
human_kinase2 = ['O14746',  # TERT
                 'FAKEPRCRAVRSLL', '12', 'S', '10224060', 'PKB_group',
                 'LTP', 'Homo sapiens', '2004-12-31 00:00:00+01']

raw_data = [columns, non_human_no_kinase, human_no_kinase, human_kinase1,
            human_kinase2]


def test_json_processing():
    test_json = _get_json_from_entry_rows(iter(raw_data))
    assert len(test_json) == len(raw_data) - 1, len(test_json)
    assert set(columns) == set(test_json[0].keys())


def test_keep_empty():
    pep = PhosphoElmProcessor(
        phosphoelm_data=_get_json_from_entry_rows(iter(raw_data))
    )
    pep.process_phosphorylations(skip_empty=False)
    stmts = pep.statements
    assert len(stmts) == 3, len(stmts)
    assert all(isinstance(st, Phosphorylation) for st in stmts)
    assert all(st.evidence[0].source_api == 'phosphoelm' for st in stmts)
    assert all(len(st.evidence[0].annotations['sequence']) > 0
               for st in stmts)


def test_skip_empty():
    pep = PhosphoElmProcessor(
        phosphoelm_data=_get_json_from_entry_rows(iter(raw_data))
    )
    pep.process_phosphorylations(skip_empty=True)
    stmts = pep.statements
    assert len(stmts) == 2, len(stmts)


def test_special_cases():
    # See http://phospho.elm.eu.org/kinases.html for list of kinases
    ag = _agent_from_str('Aurora A')
    assert ag.db_refs.get('HGNC') == '11393'

    ag = _agent_from_str('CCDPK')  # is non-human
    assert ag is None, ag

    ag = _agent_from_str('MAP2K_group')
    assert ag.db_refs.get('FPLX') == 'MAP2K'

    ag = _agent_from_str('PDHK4')
    assert ag.db_refs.get('HGNC') == '8812'

    # PDKC is probably a typo in the kinase table at
    # http://phospho.elm.eu.org/kinases.html but it is unclear what was
    # meant by the name from the source material: possibly PKC or SDK1.
    ag = _agent_from_str('PDKC')
    assert ag is None

    ag = _agent_from_str('PKA_alpha')
    assert ag.db_refs.get('HGNC') == '9380'

    ag = _agent_from_str('PKC_zeta')
    assert ag.db_refs.get('HGNC') == '9412'

    ag = _agent_from_str('RSK-5')
    assert ag.db_refs.get('HGNC') == '10434'

    ag = _agent_from_str('Titin kinase')
    assert ag.db_refs.get('HGNC') == '12403'

