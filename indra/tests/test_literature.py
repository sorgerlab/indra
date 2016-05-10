from indra.literature import id_lookup, get_full_text

def test_get_full_text_pmc():
    txt, txt_format = get_full_text('PMC4322985', 'pmcid')
    assert(txt_format == 'nxml')
    assert(len(txt) > 300000)

def test_get_full_text_doi():
    txt, txt_format = get_full_text('10.18632/oncotarget.2555', 'doi')
    assert(txt_format == 'nxml')
    assert(len(txt) > 300000)

def test_get_full_text_pubmed_abstract():
    # DOI lookup in CrossRef fails for this one because of page mismatch
    txt, txt_format = get_full_text('27075779', 'pmid')
    assert(txt_format == 'abstract')
    assert(len(txt) > 800)

def test_id_lookup():
    res = id_lookup('17513615', 'pmid')
    assert res['doi'] == '10.1158/1535-7163.mct-06-0807'

def test_id_lookup_no_pmid():
    res = id_lookup('10.1083/jcb.1974if', 'doi')
    assert res['pmcid'] == 'PMC3352949'
    res = id_lookup('PMC3352949', 'pmcid')
    assert res['doi'] == '10.1083/jcb.1974if'

