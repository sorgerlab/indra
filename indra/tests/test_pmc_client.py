from indra.literature import pmc_client

def test_get_xml():
    pmc_id = '4322985'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is not None)

def test_get_xml_PMC():
    pmc_id = 'PMC4322985'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is not None)

def test_get_xml_invalid():
    pmc_id = '9999999'
    xml_str = pmc_client.get_xml(pmc_id)
    assert(xml_str is None)
