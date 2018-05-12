from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.sources.medscan.processor import *
from indra.sources.medscan.processor import _urn_to_db_refs
from indra.sources.medscan.api import *
from indra.sources.medscan.processor import ProteinSiteInfo

# Path to the Medscan test/dummy data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this, 'medscan_tests_data')


def test_urn_to_db_refs():
    # Test resolving a subset of the urns that are used for grounding
    # The processor adds the TEXT property to the generated db_refs after
    # calling urn_to_db_refs

    # agi-cas
    urn1 = 'urn:agi-cas:89-73-6'
    db_refs_1, _ = _urn_to_db_refs(urn1)
    assert(db_refs_1 == {'CHEBI': 'CHEBI:45615'})

    # agi-llid
    urn2 = 'urn:agi-llid:9451'
    db_refs_2, hgnc_name = _urn_to_db_refs(urn2)
    assert(db_refs_2 == {'HGNC': '3255', 'UP': 'Q9NZJ5'})
    assert(hgnc_name == 'EIF2AK3')

    # agi-ncimorgan
    urn3 = 'urn:agi-ncimorgan:C0012144'
    db_refs_3, _ = _urn_to_db_refs(urn3)
    assert(db_refs_3 == {'MESH': 'C0012144'})

    # agi-nicmcelltype
    urn4 = 'urn:agi-ncimcelltype:C0242633'
    db_refs_4, _ = _urn_to_db_refs(urn4)
    assert(db_refs_4 == {'MESH': 'C0242633'})

    # agi-meshdist
    urn5 = 'urn:agi-meshdis:Paramyotonia%20Congenita'
    db_refs_5, _ = _urn_to_db_refs(urn5)
    assert(db_refs_5 == {'MESHDIS': 'Paramyotonia%20Congenita'})

    # agi-gocomplex
    urn6 = 'urn:agi-gocomplex:0005610'
    db_refs_6, _ = _urn_to_db_refs(urn6)
    assert(db_refs_6 == {'GO': 'GO:0005610', 'FPLX': 'Laminin_332'})

    # agi-go
    urn7 = 'urn:agi-go:0001515'
    db_refs_7, _ = _urn_to_db_refs(urn7)
    assert(db_refs_7 == {'GO': 'GO:0001515'})

    # agi-ncimtissue
    urn8 = 'urn:agi-ncimtissue:C0007807'
    db_refs_8, _ = _urn_to_db_refs(urn8)
    assert(db_refs_8 == {'MESH': 'C0007807'})

    # Do we ground to Famplex when there is a correspondence between a GO
    # id and a Famplex id?
    urn9 = 'urn:agi-go:0000776'
    db_refs_9, _ = _urn_to_db_refs(urn9)
    assert(db_refs_9 == {'GO': 'GO:0000776', 'FPLX': 'Kinetochore'})

    # Do we ground to Famplex when there is a correspondence between a MESH
    # id and a Famplex id?
    urn10 = 'urn:agi-ncimcelltype:D000199'
    db_refs_10, _ = _urn_to_db_refs(urn10)
    assert(db_refs_10 == {'MESH': 'D000199', 'FPLX': 'Actin'})

    # If the urn corresponds to an eccode, do we ground to famplex if that
    # eccode is in the Famplex equivalences table?
    urn11 = 'urn:agi-enz:1.1.1.1'
    db_refs_11, _ = _urn_to_db_refs(urn11)
    assert(db_refs_11 == {'FPLX': 'ADH'})

    # Do we check the Famplex equivalences table to see if a raw Medscan URN
    # maps to a Famplex ID?
    urn11 = 'urn:agi-aopfc:0000105'
    db_refs_11, _ = _urn_to_db_refs(urn11)
    assert(db_refs_11 == {'FPLX': 'GATA'})


def test_agent_from_entity():
    mp = MedscanProcessor()

    # Test entity
    entity = MedscanEntity(name='kinesin-I',
                           urn='urn:agi-gocomplex:0016938',
                           type=None, properties={})

    # Test relation
    tagged_sentence = '{ID{321=BRAF} is a protein, not a type of car.'
    relation = MedscanRelation(uri=None,
                               sec=None,
                               entities={'123': entity},
                               tagged_sentence=tagged_sentence,
                               subj=None,
                               verb=None,
                               obj=None,
                               svo_type=None)

    # Test for when an entity is in the grounded entities list
    agent1 = mp.agent_from_entity(relation, 'ID{123}')
    assert(agent1.db_refs == {'TEXT': 'kinesin-I', 'GO': 'GO:0016938'})

    # Test for when an entity is in the tagged sentence but not the entity list
    agent2 = mp.agent_from_entity(relation, 'ID{321}')
    assert(agent2.db_refs == {'TEXT': 'BRAF'})  # No grounding

    # Test for when an entity is neither tagged in the sentence nor in the
    # grounded entities list
    agent3 = mp.agent_from_entity(relation, 'ID{444}')
    assert(agent3 is None)


def test_expressioncontrol_positive():
    fname = os.path.join(data_folder, 'test_ExpressionControl_positive.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 2)

    s0 = statements[0]

    assert(s0.subj.db_refs == {'TEXT': 'hypoxia'})
    assert(s0.obj.db_refs == {'HGNC': '3415', 'TEXT': 'erythropoietin',
                              'UP': 'P01588'})


def test_evidence():
    # Test that evidence object is created correctly
    fname = os.path.join(data_folder, 'test_ExpressionControl_positive.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 2)
    s0 = statements[0]

    assert(len(s0.evidence) == 1)
    assert(isinstance(s0, IncreaseAmount))
    assert(s0.evidence[0].source_api == 'medscan')
    assert(s0.evidence[0].source_id == 'info:pmid/23455322')
    assert(s0.evidence[0].pmid == '23455322')
    assert(s0.evidence[0].text == 'Finally, we show that parp-1(-/-) mice' +
                                  ' display a significant reduction in the' +
                                  ' circulating hypoxia-induced ' +
                                  'erythropoietin levels, number of ' +
                                  'red cells and hemoglobin concentration. ')


def test_molsynthesis_positive():
    fname = os.path.join(data_folder, 'test_MolSynthesis-positive.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, IncreaseAmount))

    assert(s0.subj.db_refs == {'HGNC': '19260', 'TEXT': 'BLT2',
                               'UP': 'Q9NPC1'})
    assert(s0.obj.db_refs == {'TEXT': 'reactive oxygen species'})


def test_expressioncontrol_negative():
    fname = os.path.join(data_folder, 'test_ExpressionControl_negative.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, DecreaseAmount))
    assert(s0.subj.db_refs == {'CHEBI': 'CHEBI:6700', 'TEXT': 'matrine'})
    assert(s0.obj.db_refs == {'HGNC': '6364',
                              'TEXT': 'PSA and androgen receptor',
                              'UP': 'P07288'})


def test_molsynthesis_negative():
    fname = os.path.join(data_folder, 'test_MolSynthesis-negative.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, DecreaseAmount))

    assert(s0.subj.db_refs == {'HGNC': '9070', 'TEXT': 'pleckstrin',
                               'UP': 'P08567'})
    assert(s0.obj.db_refs == {'CHEBI': 'CHEBI:16595', 'TEXT': 'Ins(1,4,5)P3'})


def test_binding():
    fname = os.path.join(data_folder, 'test_Binding.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, Complex))
    members = s0.members
    assert(len(members) == 2)
    m0 = members[0]
    m1 = members[1]

    assert(m0.db_refs == {'HGNC': '7664', 'TEXT': 'Both Nck and Grb4',
                          'UP': 'P16333'})
    assert(m1.db_refs == {'HGNC': '9406', 'TEXT': 'PRK2', 'UP': 'Q16513'})


def test_phosphorylate():
    fname = os.path.join(data_folder, 'test_Phosphorylate.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, Phosphorylation))

    assert(s0.enz.db_refs == {'GO': 'GO:0005610', 'FPLX': 'Laminin_332',
                              'TEXT': 'IKK alpha'})
    assert(s0.enz.name == 'Laminin_332')  # agent name is FPLX when available
    assert(s0.sub.db_refs == {'HGNC': '6120', 'TEXT': 'IRF-5', 'UP': 'Q13568'})


def test_activation():
    fname = os.path.join(data_folder, 'test_Activation.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(type(s0) == Activation)
    assert(s0.subj.name == 'Laminin_332')
    assert(s0.obj.name == 'IRF-5 dimers')


def test_inhibition():
    fname = os.path.join(data_folder, 'test_Inhibition.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(type(s0) == Inhibition)
    assert(repr(s0.subj) == 'DNMT3A(mods: (methylation, R, 882))')
    assert(repr(s0.obj) == 'cell differentiation()')


def test_dephosphorylate():
    fname = os.path.join(data_folder, 'test_Dephosphorylate.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, Dephosphorylation))

    assert(s0.enz.db_refs == {'HGNC': '30579', 'TEXT': 'Slingshot-1 (SSH1',
                              'UP': 'Q8WYL5'})
    assert(s0.sub.db_refs == {'HGNC': '1874', 'TEXT': 'cofilin',
                              'UP': 'P23528'})


def test_protein_mutation():
    fname = os.path.join(data_folder, 'test_Protein_Mutation.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, Complex))

    members = s0.members
    assert(len(members) == 2)
    m0 = members[0]
    m1 = members[1]

    assert(m0.db_refs == {'HGNC': '7910', 'UP': 'P06748', 'TEXT': 'NPM1'})
    assert(len(m0.mods) == 0)
    assert(len(m0.mutations) == 0)

    assert(m1.db_refs == {'HGNC': '25994', 'UP': 'Q08J23', 'TEXT': 'NSUN2'})
    assert(len(m1.mods) == 0)
    assert(len(m1.mutations) == 1)
    mut = m1.mutations[0]
    assert(mut.position == '139')
    assert(mut.residue_from == 'S')
    assert(mut.residue_to == 'A')


def test_protein_methsite():
    fname = os.path.join(data_folder, 'test_Protein_MethSite.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, Complex))

    members = s0.members
    assert(len(members) == 2)
    m0 = members[0]
    m1 = members[1]

    assert(m0.db_refs == {'HGNC': '2978', 'UP': 'Q9Y6K1', 'TEXT': 'DNMT3A'})
    assert(len(m0.mutations) == 0)
    assert(len(m0.mods) == 1)
    mod = m0.mods[0]
    assert(mod.mod_type == 'methylation')
    assert(mod.residue == 'R')
    assert(mod.position == '882')

    assert(len(m1.mutations) == 0)
    assert(len(m1.mods) == 0)


def test_protein_phosphosite():
    fname = os.path.join(data_folder, 'test_Protein_PhosphoSite.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(isinstance(s0, DecreaseAmount))

    subj = s0.subj
    assert(subj.db_refs == {'HGNC': '115', 'UP': 'P53396', 'TEXT': 'ACLY'})
    assert(len(subj.mutations) == 0)
    assert(len(subj.mods) == 1)
    mod = subj.mods[0]
    assert(mod.residue == 'S')
    assert(mod.position == '455')
    assert(mod.mod_type == 'phosphorylation')

    obj = s0.obj
    assert(obj.db_refs == {'CHEBI': 'CHEBI:15351', 'TEXT': 'acetyl-CoA'})
    assert(len(obj.mutations) == 0)
    assert(len(obj.mods) == 0)


def test_handle_duplicates():
    # Does the processor detect duplicate SVOs within the same sentence?
    fname = os.path.join(data_folder, 'test_duplicate_SVO.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)


def test_modification_site():
    # Can we detect the modification site and residue in a modification
    # event?
    fname = os.path.join(data_folder, 'test_modification_site.csxml')
    mp = process_file(fname, None)

    statements = mp.statements
    assert(len(statements) == 1)

    s0 = statements[0]
    assert(s0.residue == 'K')
    assert(s0.position == '8')


def test_site_text_parser():
    si = ProteinSiteInfo('S10 and S20 residues', None)
    sites = si.get_sites()
    assert(len(sites) == 2)
    assert(sites[0].residue == 'S')
    assert(sites[0].position == '10')
    assert(sites[1].residue == 'S')
    assert(sites[1].position == '20')
