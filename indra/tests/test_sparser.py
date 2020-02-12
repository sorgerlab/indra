import json
from indra.sources import sparser
from indra.sources.sparser.processor import fix_agent
from indra.statements import Agent, Phosphorylation, Complex


def test_fix_agent_be_name():
    a = Agent('XXX', db_refs={'FPLX': 'CDK'})
    fix_agent(a)
    assert a.name == 'CDK'


def test_fix_agent_hgnc_only():
    a = Agent('XXX', db_refs={'HGNC': '7199'})
    fix_agent(a)
    assert a.name == 'MOS'
    assert a.db_refs.get('UP') == 'P00540'


def test_fix_agent_fa_only():
    a = Agent('XXX', db_refs={'FA': '00815'})
    fix_agent(a)
    assert a.name == 'Cyclin'
    assert a.db_refs.get('FPLX') == 'Cyclin'
    assert a.db_refs.get('NXPFA') == '00815'
    assert 'FA' not in a.db_refs


def test_fix_agent_ncit_only():
    a = Agent('XXX', db_refs={'NCIT': 'C25785'})
    fix_agent(a)
    assert a.name == 'KRAS'
    assert a.db_refs.get('HGNC') == '6407'
    assert a.db_refs.get('UP') == 'P01116'


def test_fix_agent_ncit_only():
    a = Agent('XXX', db_refs={'NCIT': 'C129655'})
    fix_agent(a)
    assert a.name == 'TUBB'
    assert a.db_refs.get('FPLX') == 'TUBB'


def test_fix_agent_pcid():
    a = Agent('XXX', db_refs={'PCID': '123'})
    fix_agent(a)
    assert 'PCID' not in a.db_refs
    assert a.db_refs['PUBCHEM'] == '123'


# ############################
# JSON processing tests
# ############################

def test_process_json_str():
    sp = sparser.process_json_dict(json.loads(json_str1))
    assert sp is not None
    assert len(sp.statements) == 1
    assert isinstance(sp.statements[0], Phosphorylation)
    sp.set_statements_pmid('1234567')
    assert sp.statements[0].evidence[0].pmid == '1234567'
    assert sp.json_stmts[0]['evidence'][0]['pmid'] == '1234567'


def test_process_json_str_with_bad_agents():
    sp = sparser.process_json_dict(json.loads(json_str2))
    assert sp is not None
    assert len(sp.statements) == 2, len(sp.statements)
    types = {type(s) for s in sp.statements}
    assert types == {Complex, Phosphorylation}, types
    assert all(len(s.agent_list()) == 2 for s in sp.statements)


def test_process_json_str_with_missing_agent():
    """This makes sure an error isn't raised in this case."""
    sp = sparser.process_json_dict(json.loads(json_str3))
    assert sp is not None
    assert len(sp.statements) == 0, len(sp.statements)


json_str1 = '''
[
 {
  "type": "Phosphorylation",
  "evidence": [
  {
    "source_api": "sparser",
    "text": "MEK phosphorylates ERK",
    "pmid": "PMC_3500"}],
  "sub": {
    "name": "ERK",
    "db_refs": {
      "NCIT": "C26360",
      "TEXT": "ERK"},
    "TEXT": "ERK"},
  "enz": {
    "name": "MEK",
    "db_refs": {
      "FPLX": "MEK",
      "TEXT": "MEK"},
    "TEXT": "MEK"}
 }
]'''

json_str2 = '''
[
  {
    "type": "Phosphorylation",
    "evidence": [
      {
          "source_api": "sparser",
          "text": "MEK phosphorylates ERK",
          "pmid": "PMC_3500"}],
    "sub": "ERK",
    "enz": {
      "name": "MEK",
      "db_refs": {
          "FPLX": "MEK",
          "TEXT": "MEK"},
      "TEXT": "MEK"}
  },
  {
    "type": "Complex",
    "members": [
      "MEK",
      {
        "name": "ERK",
        "db_refs": {"FPLX": "ERK",
                    "TEXT": "ERK"}
      }
    ],
    "belief": 1,
    "id": "3eedc7a9-fbbd-4e2e-b227-07d96f4bcff5"
  }
]'''


json_str3 = '''
[
    {
      "type": "Inhibition",
      "obj_activity": "activity",
      "evidence": [
        {
          "text": "The in vivo and in vitro studies suggested that NR enzyme is inhibited by NO in a mediated process that requires the cell integrity.",
          "source_api": "sparser",
          "pmid": "PMC10191200"
        }
      ],
      "obj": {
        "db_refs": {
          "UP": "P22945"
        },
        "name": "NIA_EMENI",
        "TEXT": "NR"
      }
    }
]
'''
