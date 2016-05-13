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
    """Look up a paper that has a PMCID and DOI but not PMID."""
    res = id_lookup('10.1083/jcb.1974if', 'doi')
    assert res['pmcid'] == 'PMC3352949'
    res = id_lookup('PMC3352949', 'pmcid')
    assert res['doi'] == '10.1083/jcb.1974if'

def test_cr_fulltext_elsevier():
    """Test the ability to obtain publications from Elsevier using the
    CrossRef Clickthrough API. Note: requires a cr_clickthrough_key file
    and an account in which the Elsevier ClickThrough License has been
    accepted."""
    # Note that this article is open access so it doesn't test the use of the
    # CrossRef clickthrough key
    (content, type) = get_full_text('23337888', 'pmid',
                                    preferred_content_type='text/xml')
    assert type == 'text/xml'
    assert len(content) == 117093
    (content, type) = get_full_text('23337888', 'pmid',
                                    preferred_content_type='text/plain')
    assert type == 'text/plain'
    assert len(content) == 53961
    # Try a smattering of other papers
    (content, type) = get_full_text('19909739', 'pmid',
                                    preferred_content_type='text/xml')
    assert type == 'text/xml'

def test_cr_fulltext_wiley():
    """Test the ability to obtain publications from Wiley using the
    CrossRef Clickthrough API. Note: requires a cr_clickthrough_key file
    and an account in which the Wiley ClickThrough License has been
    accepted."""
    (content, type) = get_full_text('20840664', 'pmid')
    assert type == 'application/pdf'
    (content, type) = get_full_text('16619251', 'pmid')
    assert type == 'application/pdf'
    # Returns 403
    (content, type) = get_full_text('12811820', 'pmid')
    assert content is None
    (content, type) = get_full_text('20803551', 'pmid')
    assert content is None

def test_other_fulltexts_with_link():
    """Test the ability to obtain publications from other publishers that have
    a full text link and a text-mining compatible license (e.g., Hindawi, which
    uses a Creative Commons License)."""
    pass

def test_fulltext_asbmb():
    (content, type) = get_full_text('14761976', 'pmid')

if __name__ == '__main__':
    test_fulltext_asbmb()

