import os
from indra.sources import biopax
import indra.sources.biopax.processor as bpc
from indra.databases import uniprot_client
from indra.util import unicode_strs
from nose.plugins.attrib import attr

model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'biopax_test.owl')


bp = biopax.process_owl(model_path)
uri_prefix = 'http://purl.org/pc2/7/'


def test_listify():
    assert bpc._listify(1) == [1]
    assert bpc._listify([1,2] == [1,2])
    assert bpc._listify([1] == [1])


def test_uniprot_id_pe():
    bpe = bp.model.objects['http://identifiers.org/reactome/REACT_117886.3']
    ids = bp._get_uniprot_id(bpe)
    assert 'Q15303' == ids


def test_name_uniprot_id_no_hgnc_pe():
    bpe = bp.model.objects['http://identifiers.org/reactome/REACT_27053.2']
    name = bp._get_element_name(bpe)
    assert name == 'IGHV3-13'


def test_uniprot_id_er():
    bpe = bp.model.objects['http://identifiers.org/uniprot/Q15303']
    ids = bp._get_uniprot_id(bpe)
    assert 'Q15303' == ids


def test_get_hgnc_id():
    bpe = bp.model.objects['http://identifiers.org/uniprot/Q15303']
    hgnc_id = bp._get_hgnc_id(bpe)
    assert hgnc_id == '3432'


def test_get_hgnc_name():
    hgnc_name = bp._get_hgnc_name('3432')
    assert hgnc_name == 'ERBB4'


def test_get_mod_feature():
    bpe = bp.model.objects[
            'ModificationFeature_bd27a53570fb9a5094bb5929bd973217']
    mc = bpc.BiopaxProcessor.mod_condition_from_mod_feature(bpe)
    assert mc.mod_type == 'phosphorylation'
    assert mc.residue == 'T'
    assert mc.position == '274'


def test_get_entity_mods():
    bpe = bp.model.objects['Protein_7aeb1631f64e49491b7a0303aaaec536']
    mods = bpc.BiopaxProcessor._get_entity_mods(bpe)
    assert len(mods) == 5
    mod_types = set([m.mod_type for m in mods])
    assert mod_types == set(['phosphorylation'])
    residues = set([m.residue for m in mods])
    assert residues == set(['Y'])
    mod_pos = set([m.position for m in mods])
    assert mod_pos == set(['1035', '1056', '1128', '1188', '1242'])


@attr('webservice', 'slow')
def test_pathsfromto():
    bp = biopax.process_pc_pathsfromto(['MAP2K1'], ['MAPK1'])
    assert_pmids(bp.statements)
    assert_source_sub_id(bp.statements)
    assert unicode_strs(bp.statements)
    num_unique = len({s.get_hash(shallow=False) for s in bp.statements})
    assert len(bp.statements) == num_unique


def test_all_uniprot_ids():
    for obj in bp.model.objects.values():
        if bpc._is_protein(obj):
            uniprot_id = bp._get_uniprot_id(obj)
            if uniprot_id is not None:
                assert not uniprot_client.is_secondary(uniprot_id)
                assert unicode_strs(uniprot_id)


def test_all_hgnc_ids():
    for obj in bp.model.objects.values():
        if bpc._is_protein(obj):
            hgnc_id = bp._get_hgnc_id(obj)
            if hgnc_id is not None:
                assert unicode_strs(hgnc_id)


def test_all_protein_db_refs():
    unmapped_uniprot_ids = []
    for obj in bp.model.objects.values():
        if bpc._is_protein(obj):
            db_refs = bpc.BiopaxProcessor._get_db_refs(obj)
            uniprot_id = db_refs.get('UP')
            hgnc_id = db_refs.get('HGNC')
            if uniprot_id:
                if uniprot_client.is_human(uniprot_id):
                    if not hgnc_id:
                        unmapped_uniprot_ids.append(uniprot_id)
    unmapped_uniprot_ids = sorted(list(set(unmapped_uniprot_ids)))
    # The number of unmapped entries should not increase
    # so we check for an upper limit here
    assert len(unmapped_uniprot_ids) < 95


def assert_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                assert ev.pmid.isdigit()


def assert_source_sub_id(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            assert 'source_sub_id' in ev.annotations
            assert ev.annotations['source_sub_id']
