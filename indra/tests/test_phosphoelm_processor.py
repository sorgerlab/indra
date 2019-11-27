import csv
from indra.util.aws import get_s3_client
from indra.statements import Phosphorylation
from indra.sources.phosphoelm.api import PhosphoElmProcessor, \
    _get_json_from_entry_rows, s3_bucket, ppelm_s3_key
from nose.plugins.attrib import attr

columns = ['acc', 'sequence', 'position', 'code', 'pmids', 'kinases',
           'source', 'species', 'entry_date']
non_human_no_kinase = ['O08539',
                       'FAKEGSKG', '6', 'S', '17114649', '', 'HTP',
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


@attr('nonpublic')
def test_download_from_s3():
    s3 = get_s3_client(False)
    s3_obj = s3.get_object(Bucket=s3_bucket, Key=ppelm_s3_key)
    csv_reader = csv.reader(
        s3_obj['Body'].read().decode('utf8').splitlines(True),
        delimiter='\t'
    )
    ppelm_json = _get_json_from_entry_rows(csv_reader)
    assert ppelm_json


def test_json_processing():
    test_json = _get_json_from_entry_rows(iter(raw_data))
    assert len(test_json) == len(raw_data) - 1, len(test_json)
    assert set(columns) == set(test_json[0].keys())


def test_keep_empty():
    pep = PhosphoElmProcessor(
        phosphoelm_data=_get_json_from_entry_rows(iter(raw_data))
    )
    pep.process_phosphorylations(True)
    stmts = pep.statements
    assert len(stmts) == 3
    assert all(isinstance(st, Phosphorylation) for st in stmts)
    assert all(st.evidence[0].source_api == 'phospho.ELM' for st in stmts)
    assert all(len(st.evidence[0].annotations['sequence']) > 0
               for st in stmts)


def test_not_empty():
    pep = PhosphoElmProcessor(
        phosphoelm_data=_get_json_from_entry_rows(iter(raw_data))
    )
    pep.process_phosphorylations()
    stmts = pep.statements
    assert len(stmts) == 2
