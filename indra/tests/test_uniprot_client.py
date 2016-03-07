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
    assert(uniprot_client.get_hgnc_name(g) == 'MAPK3')

def test_get_family_members():
    members = uniprot_client.get_family_members('RAF')
    assert('ARAF' in members)
    assert('BRAF' in members)
    assert('RAF1' in members)

def test_get_hgnc_name_human():
    g = uniprot_client.query_protein('P00533')
    hgnc_name = uniprot_client.get_hgnc_name(g)
    assert(hgnc_name == 'EGFR')

def test_get_hgnc_name_nonhuman():
    g = uniprot_client.query_protein('P31938')
    hgnc_name = uniprot_client.get_hgnc_name(g)
    assert(hgnc_name is None)

def test_get_gene_name_human():
    g = uniprot_client.query_protein('P00533')
    gene_name = uniprot_client.get_gene_name(g)
    assert(gene_name == 'EGFR')

def test_get_gene_name_nonhuman():
    g = uniprot_client.query_protein('P31938')
    gene_name = uniprot_client.get_gene_name(g)
    assert(gene_name == 'Map2k1')

def test_get_sequence():
    g = uniprot_client.query_protein('P00533')
    seq = uniprot_client.get_sequence(g)
    assert(len(seq) > 1000)

def test_get_modifications():
    g = uniprot_client.query_protein('P27361')
    mods = uniprot_client.get_modifications(g)
    assert(('Phosphothreonine', 202) in mods)
    assert(('Phosphotyrosine', 204) in mods)

def test_verify_location():
    g = uniprot_client.query_protein('P27361')
    assert(uniprot_client.verify_location(g, 'T', 202)) 
    assert(not uniprot_client.verify_location(g, 'S', 202))
    assert(not uniprot_client.verify_location(g, 'T', -1))
    assert(not uniprot_client.verify_location(g, 'T', 10000))
