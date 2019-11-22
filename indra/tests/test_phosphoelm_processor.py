from indra.statements import Phosphorylation
from indra.sources.phosphoELM.api import PhosphoELMPRocessor, \
    _get_json_from_entry_rows

columns = ['acc', 'sequence', 'position', 'code', 'pmids', 'kinases', 'source',
           'species', 'entry_date']
non_human_no_kinase = ['O08539',
                       'MAEMGSKG', '6', 'S', '17114649', '', 'HTP',
                       'Mus musculus', '2005-03-14 12:16:11.108314+01']
human_no_kinase = ['O14543',
                   'MVTHSKFPAAGMSRPLDTSLRLKTFSSKSEYQL', '31', 'Y',
                   '12783885', '', 'LTP', 'Homo sapiens',
                   '2006-10-17 12:06:48.271076+02']
human_kinase1 = ['O14543',
                 'FAKESEQUENCESRPLDTSLRLKTFSSKSEYQL',
                 '31', 'Y', '12783885', 'Lck', 'LTP', 'Homo sapiens',
                 '2006-10-17 12:06:48.16767+02']
human_kinase2 = ['O14746', 'FAKEPRCRAVRSLL', '12', 'S', '10224060',
                 'PKB_group', 'LTP', 'Homo sapiens',
                 '2004-12-31 00:00:00+01']

raw_data = [columns, non_human_no_kinase, human_no_kinase, human_kinase1,
            human_kinase2]


def test_json_processing():
    test_json = _get_json_from_entry_rows(iter(raw_data))
    assert len(test_json) == len(raw_data) - 1, len(test_json)
    assert set(columns) == set(test_json[0].keys())


def test_keep_empty():
    stmts = PhosphoELMPRocessor(
        file_dump_json=_get_json_from_entry_rows(iter(raw_data)),
        keep_empty=True
    ).statements
    assert len(stmts) == 3
    assert all(isinstance(st, Phosphorylation) for st in stmts)
    assert all(st.evidence[0].source_api == 'phospho.ELM' for st in stmts)
    assert all(len(st.evidence[0].annotations['sequence']) > 0
               for st in stmts)


def test_not_empty():
    stmts = PhosphoELMPRocessor(
        file_dump_json=_get_json_from_entry_rows(iter(raw_data))
    ).statements
    assert len(stmts) == 2


