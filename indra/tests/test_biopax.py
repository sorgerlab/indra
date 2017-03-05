from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
from indra.java_vm import autoclass, cast
from indra import biopax
import indra.biopax.processor as bpc
from indra.databases import uniprot_client
from indra.util import unicode_strs
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies

model_path = os.path.dirname(os.path.abspath(__file__)) +\
             '/../../data/biopax_test.owl'

bp = biopax.process_owl(model_path)
uri_prefix = 'http://purl.org/pc2/7/'

def test_paxtools_autoclass():
    autoclass('org.biopax.paxtools.impl.level3.ProteinImpl')

def test_biopaxpattern_autoclass():
    autoclass('org.biopax.paxtools.pattern.PatternBox')

def test_cpath_autoclass():
    autoclass('cpath.client.CPathClient')

def test_listify():
    assert(bpc._listify(1) == [1])
    assert(bpc._listify([1,2] == [1,2]))
    assert(bpc._listify([1] == [1]))

def test_list_listify():
    assert(bpc._list_listify([1]) == [[1]])
    assert(bpc._list_listify([1,2]) == [[1],[2]])
    assert(bpc._list_listify([1, [1,2]]) == [[1], [1,2]])

def test_get_combinations():
    combs = [c for c in bpc._get_combinations([1, 2])]
    assert(combs == [(1,2)])
    combs = [c for c in bpc._get_combinations([1, [3,4]])]
    assert(combs == [(1,3), (1,4)])

def test_has_members_er():
    bpe = bp.model.getByID(uri_prefix +\
                     'ProteinReference_971cec47bcd850e2b7d602f0416edacf')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    assert(bpc._has_members(bpe))

    bpe = bp.model.getByID('http://identifiers.org/uniprot/P56159')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    assert(not bpc._has_members(bpe))

def test_has_members_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc._bp('Protein'), bpe)
    assert(bpc._has_members(bpe))

def test_has_members_pe2():
    bpe = bp.model.getByID(uri_prefix + 'Protein_7d526475fd43d0a07ca1a596fe81aae0')
    bpe = cast(bpc._bp('Protein'), bpe)
    assert(not bpc._has_members(bpe))

def test_is_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc._bp('Protein'), bpe)
    assert(bpc._is_entity(bpe))

def test_is_pe2():
    bpe = bp.model.getByID(uri_prefix +\
                     'ProteinReference_971cec47bcd850e2b7d602f0416edacf')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    assert(not bpc._is_entity(bpe))

def test_is_er():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc._bp('Protein'), bpe)
    assert(not bpc._is_reference(bpe))

def test_is_er2():
    bpe = bp.model.getByID(uri_prefix +\
                     'ProteinReference_971cec47bcd850e2b7d602f0416edacf')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    assert(bpc._is_reference(bpe))

def test_is_mod():
    bpe = bp.model.getByID(uri_prefix +\
                    'ModificationFeature_59c99eae672d2a11e971a93c7848d5c6')
    bpe = cast(bpc._bp('ModificationFeature'), bpe)
    assert(bpc._is_modification(bpe))

def test_is_mod2():
    bpe = bp.model.getByID(uri_prefix +\
                    'FragmentFeature_806ae27c773eb2d9138269552899c242')
    bpe = cast(bpc._bp('FragmentFeature'), bpe)
    assert(not bpc._is_modification(bpe))

def test_is_complex():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_24213.2')
    bpe = cast(bpc._bp('Complex'), bpe)
    assert(bpc._is_complex(bpe))

def test_is_complex2():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117345.2')
    bpe = cast(bpc._bp('Protein'), bpe)
    assert(not bpc._is_complex(bpe))

def test_uniprot_id_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_117886.3')
    bpe = cast(bpc._bp('Protein'), bpe)
    ids = bp._get_uniprot_id(bpe)
    assert('Q15303' == ids)

def test_name_uniprot_id_no_hgnc_pe():
    bpe = bp.model.getByID('http://identifiers.org/reactome/REACT_27053.2')
    bpe = cast(bpc._bp('Protein'), bpe)
    name = bp._get_element_name(bpe)
    assert(name == 'IGHV3-13')

def test_uniprot_id_er():
    bpe = bp.model.getByID('http://identifiers.org/uniprot/Q15303')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    ids = bp._get_uniprot_id(bpe)
    assert('Q15303' == ids)

def test_get_hgnc_id():
    bpe = bp.model.getByID('http://identifiers.org/uniprot/Q15303')
    bpe = cast(bpc._bp('ProteinReference'), bpe)
    hgnc_id = bp._get_hgnc_id(bpe)
    assert(hgnc_id == '3432')

def test_get_hgnc_name():
    hgnc_name = bp._get_hgnc_name('3432')
    assert(hgnc_name == 'ERBB4')

def test_get_mod_feature():
    bpe = bp.model.getByID(uri_prefix +\
            'ModificationFeature_bd27a53570fb9a5094bb5929bd973217')
    mf = cast(bpc._bp('ModificationFeature'), bpe)
    mc = bpc.BiopaxProcessor._extract_mod_from_feature(mf)
    assert(mc[0] == 'phosphorylation')
    assert(mc[1] == 'T')
    assert(mc[2] == '274')

def test_get_entity_mods():
    bpe = bp.model.getByID(uri_prefix +\
            'Protein_7aeb1631f64e49491b7a0303aaaec536')
    protein = cast(bpc._bp('Protein'), bpe)
    mods = bpc.BiopaxProcessor._get_entity_mods(protein)
    assert(len(mods) == 5)
    mod_types = set([m[0] for m in mods])
    assert(mod_types == set(['phosphorylation']))
    residues = set([m[1] for m in mods])
    assert(residues == set(['Y']))
    mod_pos = set([m[2] for m in mods])
    assert(mod_pos == set(['1035', '1056', '1128', '1188', '1242']))

def test_pathsfromto():
    bp = biopax.process_pc_pathsfromto(['MAP2K1'], ['MAPK1'])
    bp.get_phosphorylation()
    assert_pmids(bp.statements)
    assert(unicode_strs(bp.statements))

def test_all_uniprot_ids():
    for obj in bp.model.getObjects().toArray():
        bpe = bpc._cast_biopax_element(obj)
        if bpc._is_protein(bpe):
            uniprot_id = bp._get_uniprot_id(bpe)
            if uniprot_id is not None:
                assert(not uniprot_client.is_secondary(uniprot_id))
                assert(unicode_strs(uniprot_id))

def test_all_hgnc_ids():
    for obj in bp.model.getObjects().toArray():
        bpe = bpc._cast_biopax_element(obj)
        if bpc._is_protein(bpe):
            hgnc_id = bp._get_hgnc_id(bpe)
            if hgnc_id is not None:
                assert(unicode_strs(hgnc_id))

def test_all_protein_db_refs():
    unmapped_uniprot_ids = []
    for obj in bp.model.getObjects().toArray():
        bpe = bpc._cast_biopax_element(obj)
        if bpc._is_protein(bpe):
            db_refs = bpc.BiopaxProcessor._get_db_refs(bpe)
            uniprot_id = db_refs.get('UP')
            hgnc_id = db_refs.get('HGNC')
            if uniprot_id:
                if uniprot_client.is_human(uniprot_id):
                    if not hgnc_id:
                        unmapped_uniprot_ids.append(uniprot_id)
    unmapped_uniprot_ids = sorted(list(set(unmapped_uniprot_ids)))
    # The number of unmapped entries should not increase
    # so we check for an upper limit here
    assert(len(unmapped_uniprot_ids) < 95)

def assert_pmids(stmts):
    for stmt in stmts:
        for ev in stmt.evidence:
            if ev.pmid is not None:
                assert(ev.pmid.isdigit())
