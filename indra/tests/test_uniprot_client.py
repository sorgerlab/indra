from indra.databases import uniprot_client

def test_query_protein_exists():
    g = uniprot_client.query_protein('P00533')
    assert(g is not None)

def test_query_protein_nonexist():
    g = uniprot_client.query_protein('XXXX')
    assert(g is None)

def test_query_protein_deprecated():
    g = uniprot_client.query_protein('Q8NHX1')
    assert(g is not None)
    assert(uniprot_client.get_hgnc_name('Q8NHX1') == 'MAPK3')

def test_get_family_members():
    members = uniprot_client.get_family_members('RAF')
    assert('ARAF' in members)
    assert('BRAF' in members)
    assert('RAF1' in members)

def test_get_hgnc_name_human():
    hgnc_name = uniprot_client.get_hgnc_name('P00533')
    assert(hgnc_name == 'EGFR')

def test_get_hgnc_name_nonhuman():
    hgnc_name = uniprot_client.get_hgnc_name('P31938')
    assert(hgnc_name is None)

def test_get_gene_name_human():
    gene_name = uniprot_client.get_gene_name('P00533')
    assert(gene_name == 'EGFR')

def test_get_gene_name_nonhuman():
    gene_name = uniprot_client.get_gene_name('P31938')
    assert(gene_name == 'Map2k1')

def test_get_sequence():
    seq = uniprot_client.get_sequence('P00533')
    assert(len(seq) > 1000)

def test_get_modifications():
    mods = uniprot_client.get_modifications('P27361')
    assert(('Phosphothreonine', 202) in mods)
    assert(('Phosphotyrosine', 204) in mods)

def test_verify_location():
    assert(uniprot_client.verify_location('P27361', 'T', 202)) 
    assert(not uniprot_client.verify_location('P27361', 'S', 202))
    assert(not uniprot_client.verify_location('P27361', 'T', -1))
    assert(not uniprot_client.verify_location('P27361', 'T', 10000))

def test_get_mnemonic():
    mnemonic = uniprot_client.get_mnemonic('Q02750')
    assert(mnemonic == 'MP2K1_HUMAN')
