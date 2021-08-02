import re
from indra.databases.identifiers import get_identifiers_url, \
    parse_identifiers_url, get_ns_from_identifiers,\
    get_ns_id_from_identifiers, get_identifiers_ns, namespace_embedded, \
    ensure_prefix_if_needed, identifiers_registry


def test_map_ns():
    assert get_ns_from_identifiers('go') == 'GO'
    assert get_ns_from_identifiers('uniprot') == 'UP'
    assert get_ns_from_identifiers('XXX') is None


def test_map_ns_id():
    assert get_ns_id_from_identifiers('uniprot', 'P12345') == \
        ('UP', 'P12345')
    assert get_ns_id_from_identifiers('go', 'GO:0005856') == \
        ('GO', 'GO:0005856')


def test_identifiers_ns():
    assert get_identifiers_ns('UP') == 'uniprot'
    assert get_identifiers_ns('GO') == 'go'
    assert get_identifiers_ns('XXX') is None


# For each pair of namespace and ID we get all possible URL forms.
# First URL in each tuple corresponds to the latest standard format.
ns_mapping = {
    ('UP', 'P0DP23'): ('https://identifiers.org/uniprot:P0DP23',
                       'http://identifiers.org/uniprot:P0DP23',
                       'https://identifiers.org/uniprot/P0DP23',
                       'http://identifiers.org/uniprot/P0DP23'),
    ('HGNC', '2674'): ('https://identifiers.org/hgnc:2674',
                       'http://identifiers.org/hgnc:2674',
                       'https://identifiers.org/hgnc/HGNC:2674',
                       'http://identifiers.org/hgnc/HGNC:2674'),
    ('IP', 'IPR000100'): ('https://identifiers.org/interpro:IPR000100',
                          'http://identifiers.org/interpro:IPR000100',
                          'https://identifiers.org/interpro/IPR000100',
                          'http://identifiers.org/interpro/IPR000100'),
    ('CHEBI', 'CHEBI:36927'): ('https://identifiers.org/CHEBI:36927',
                               'http://identifiers.org/CHEBI:36927',
                               'https://identifiers.org/chebi/CHEBI:36927',
                               'http://identifiers.org/chebi/CHEBI:36927'),
    ('NCIT', 'C80519'): ('https://identifiers.org/ncit:C80519',
                         'http://identifiers.org/ncit:C80519',
                         'https://identifiers.org/ncit/C80519',
                         'http://identifiers.org/ncit/C80519'),
    ('GO', 'GO:0006915'): ('https://identifiers.org/GO:0006915',
                           'http://identifiers.org/GO:0006915',
                           'https://identifiers.org/go/GO:0006915',
                           'http://identifiers.org/go/GO:0006915'),
    ('PUBCHEM', '100101'): ('https://identifiers.org/pubchem.compound:100101',
                            'http://identifiers.org/pubchem.compound:100101',
                            'https://identifiers.org/pubchem.compound/100101',
                            'http://identifiers.org/pubchem.compound/100101'),
    ('PF', 'PF01234'): ('https://identifiers.org/pfam:PF01234',
                        'http://identifiers.org/pfam:PF01234',
                        'https://identifiers.org/pfam/PF01234',
                        'http://identifiers.org/pfam/PF01234'),
    ('MIRBASEM', 'MIMAT0000001'): (
        'https://identifiers.org/mirbase.mature:MIMAT0000001',
        'http://identifiers.org/mirbase.mature:MIMAT0000001',
        'https://identifiers.org/mirbase.mature/MIMAT0000001',
        'http://identifiers.org/mirbase.mature/MIMAT0000001'),
    ('MIRBASE', 'MI0000001'): ('https://identifiers.org/mirbase:MI0000001',
                               'http://identifiers.org/mirbase:MI0000001',
                               'https://identifiers.org/mirbase/MI0000001',
                               'http://identifiers.org/mirbase/MI0000001'),
    ('MESH', 'C000100'): ('https://identifiers.org/mesh:C000100',
                          'http://identifiers.org/mesh:C000100',
                          'https://identifiers.org/mesh/C000100',
                          'http://identifiers.org/mesh/C000100'),
    ('EGID', '100010'): ('https://identifiers.org/ncbigene:100010',
                         'http://identifiers.org/ncbigene:100010',
                         'https://identifiers.org/ncbigene/100010',
                         'http://identifiers.org/ncbigene/100010'),
    ('HMDB', 'HMDB00001'): ('https://identifiers.org/hmdb:HMDB00001',
                            'http://identifiers.org/hmdb:HMDB00001',
                            'https://identifiers.org/hmdb/HMDB00001',
                            'http://identifiers.org/hmdb/HMDB00001'),
    ('FPLX', 'RAS'): ('https://identifiers.org/fplx:RAS',
                      'http://identifiers.org/fplx:RAS',
                      'https://identifiers.org/fplx/RAS',
                      'http://identifiers.org/fplx/RAS'),
    ('REFSEQ_PROT', 'NP_012345'): ('https://identifiers.org/refseq:NP_012345',
                                   'http://identifiers.org/refseq:NP_012345',
                                   'https://identifiers.org/refseq/NP_012345',
                                   'http://identifiers.org/refseq/NP_012345'),
    ('EFO', '0004859'): ('https://identifiers.org/efo:0004859',
                         'http://identifiers.org/efo:0004859',
                         'https://identifiers.org/efo/0004859',
                         'http://identifiers.org/efo/0004859'),
    ('HP', 'HP:0000118'): ('https://identifiers.org/HP:0000118',
                           'http://identifiers.org/HP:0000118',
                           'https://identifiers.org/hp/HP:0000118',
                           'http://identifiers.org/hp/HP:0000118'),
    ('DOID', 'DOID:11337'): ('https://identifiers.org/DOID:11337',
                             'http://identifiers.org/DOID:11337'),
    ('ECCODE', '1.1.1.1'): ('https://identifiers.org/ec-code:1.1.1.1',
                            'http://identifiers.org/ec-code:1.1.1.1'),
    ('CAS', '50-00-0'): ('https://identifiers.org/cas:50-00-0',
                         'http://identifiers.org/cas:50-00-0'),
    ('DRUGBANK', 'DB00001'): ('https://identifiers.org/drugbank:DB00001',
                              'http://identifiers.org/drugbank:DB00001'),
    ('TAXONOMY', '9606'): ('https://identifiers.org/taxonomy:9606',
                           'http://identifiers.org/taxonomy:9606'),
    ('BTO', 'BTO:0000146'): ('https://identifiers.org/BTO:0000146',
                             'http://identifiers.org/BTO:0000146'),
    ('CHEMBL', 'CHEMBL1229517'): (
        'https://identifiers.org/chembl.compound:CHEMBL1229517',
        'http://identifiers.org/chembl.compound:CHEMBL1229517',
        'https://identifiers.org/chembl.compound/CHEMBL1229517',
        'http://identifiers.org/chembl.compound/CHEMBL1229517'),
    ('LINCS', 'LSM-6357'): (
        'https://identifiers.org/lincs.smallmolecule:LSM-6357',
        'http://identifiers.org/lincs.smallmolecule:LSM-6357',
        'https://identifiers.org/lincs.smallmolecule/LSM-6357',
        'http://identifiers.org/lincs.smallmolecule/LSM-6357'),
    ('UPPRO', 'PRO_0000032458'): (
        'https://identifiers.org/uniprot.chain:PRO_0000032458',
        'https://identifiers.org/uniprot:P01019#PRO_0000032458',
        'http://identifiers.org/uniprot:P01019#PRO_0000032458'),
    ('NXPFA', '00880'): ('https://www.nextprot.org/term/FA-00880', ),
    ('SIGNOR', 'SIGNOR-PF15'): (
        'https://signor.uniroma2.it/relation_result.php?id=SIGNOR-PF15', ),
    ('HGNC_GROUP', '643'): (
        'https://identifiers.org/hgnc.genefamily:643', ),
    ('PR', 'PR:000000019'): (
        'https://identifiers.org/PR:000000019', ),
    ('NCBIPROTEIN', 'QHO60603'): (
        'https://identifiers.org/ncbiprotein:QHO60603', ),
    ('NONCODE', 'NONHSAT028507.2'): (
        'https://identifiers.org/noncodev4.rna:NONHSAT028507.2', ),
    ('HMS-LINCS', '10220'): ('http://lincs.hms.harvard.edu/db/sm/10220-101', ),
    ('LNCRNADB', 'URS000075C808_9606'): (
        'https://identifiers.org/rnacentral:URS000075C808_9606', ),
    ('SCHEM', None): (
        'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
        'selventa-legacy-chemicals/selventa-legacy-chemicals-20150601.belns',),
    ('SCOMP', None): (
        'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
        'selventa-named-complexes/selventa-named-complexes-20150601.belns',),
    ('SFAM', None): (
        'https://arty.scai.fraunhofer.de/artifactory/bel/namespace/'
        'selventa-protein-families/selventa-protein-families-20150601.belns',),
}


def test_chembl():
    cid = '1229517'
    assert get_identifiers_url('CHEMBL', cid) == \
        'https://identifiers.org/chembl.compound:CHEMBL%s' % cid
    assert get_identifiers_url('CHEMBL', 'CHEMBL%s' % cid) == \
        'https://identifiers.org/chembl.compound:CHEMBL%s' % cid


def test_signor():
    sid = 'SIGNOR-PF15'
    assert get_identifiers_url('SIGNOR', sid) == \
        'https://signor.uniroma2.it/relation_result.php?id=%s' % sid


def test_get_identifiers_url():
    # Get latest standard URL for a given namespace and ID
    for ns_tuple, urls in ns_mapping.items():
        url = get_identifiers_url(*ns_tuple)
        assert url == urls[0], (url, urls[0], ns_tuple)


def test_parse_identifiers_url():
    # Get correct namespace and ID from standard and outdated URL formats
    for ns_tuple, urls in ns_mapping.items():
        for url in urls:
            ns, db_id = parse_identifiers_url(url)
            assert (ns, db_id) == ns_tuple, (url, ns, db_id)


def test_namespace_embedded():
    assert namespace_embedded('CHEBI') is True
    assert namespace_embedded('GO') is True
    assert namespace_embedded('EFO') is False
    assert namespace_embedded('XXXXX') is False


def test_ensure_prefix_if_needed():
    assert ensure_prefix_if_needed('CHEBI', 'CHEBI:123') == 'CHEBI:123'
    assert ensure_prefix_if_needed('CHEBI', '123') == 'CHEBI:123'
    assert ensure_prefix_if_needed('GO', '00004') == 'GO:00004'
    assert ensure_prefix_if_needed('EFO', '1234') == '1234'
    assert ensure_prefix_if_needed('XXXX', '1234') == '1234'


def test_eccode_override():
    assert re.match(identifiers_registry['ec-code']['pattern'], '1')
