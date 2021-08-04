from indra.statements import Agent, Evidence, Translocation, Phosphorylation
from indra.tools.fix_invalidities import *
import indra.tools.assemble_corpus as ac


def test_fix_db_refs():
    assert fix_invalidities_db_refs({'CHEBI': '123'}) == \
           {'CHEBI': 'CHEBI:123'}
    assert fix_invalidities_db_refs({'CHEMBL': '123'}) == \
           {'CHEMBL': 'CHEMBL123'}
    assert fix_invalidities_db_refs({'PUBCHEM': 'CID: 123'}) == \
           {'PUBCHEM': '123'}
    assert fix_invalidities_db_refs({'ECCODE': '1.2.-.-'}) == \
           {'ECCODE': '1.2'}
    assert fix_invalidities_db_refs({'UNIPROT': 'P12345'}) == \
           {'UP': 'P12345'}
    assert fix_invalidities_db_refs({'UNIPROT': 'SL-123'}) == \
           {'UPLOC': 'SL-123'}
    assert fix_invalidities_db_refs({'MGI': 'Mgi:Abcd1'}) == \
           {}
    assert fix_invalidities_db_refs({'RGD': 'Abcd1'}) == \
           {}


def test_fix_evidence():
    ev = Evidence(pmid='123', text_refs={'pmcid': 'PMC1234'})
    fix_invalidities_evidence(ev)
    assert ev.pmid == '123'
    assert ev.text_refs == {'PMCID': 'PMC1234', 'PMID': '123'}


def test_fix_stmts():
    stmts = [Translocation(Agent('x'), to_location=None, from_location=None),
             Phosphorylation(Agent('a', db_refs={'TEXT': None,
                                                 'FPLX': 'ERK'}), Agent('b'),
                             evidence=[Evidence(text='x')])]
    stmts_out = fix_invalidities(stmts)
    assert len(stmts_out) == 1
    assert stmts_out[0].enz.db_refs == {'FPLX': 'ERK'}

    stmts_out = ac.fix_invalidities(stmts)

    assert len(stmts_out) == 1
    assert stmts_out[0].enz.db_refs == {'FPLX': 'ERK'}

    stmts_out = ac.fix_invalidities(stmts,
                                    in_place=True,
                                    print_report_before=True,
                                    print_report_after=True,
                                    prior_hash_annots=True)
    # Check the in-place effect
    assert stmts[1].enz.db_refs == {'FPLX': 'ERK'}
    assert stmts_out[0].enz.db_refs == {'FPLX': 'ERK'}

    assert stmts_out[0].evidence[0].annotations['prior_hash']