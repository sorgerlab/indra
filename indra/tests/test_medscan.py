from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import join, dirname
from nose.tools import raises

from indra.statements import *
from indra.sources.medscan.processor import *
from indra.sources.medscan.api import *

# Path to the Medscan test/dummy data folder
path_this = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(path_this, 'medscan_tests_data')

def test_urn_to_db_refs():
    # Test resolving a subset of the urns that are used for grounding
    # The processor adds the TEXT property to the generated db_refs after
    # calling urn_to_db_refs

    # agi-cas
    urn1 = 'urn:agi-cas:89-73-6'
    db_refs_1 = urn_to_db_refs(urn1)
    assert(db_refs_1 == {'CHEBI': '45615'})

    # agi-llid
    urn2 = 'urn:agi-llid:9451'
    db_refs_2 = urn_to_db_refs(urn2)
    assert(db_refs_2 == {'HGNC': '3255'})

    # agi-ncimorgan
    urn3 = 'urn:agi-ncimorgan:C0012144'
    db_refs_3 = urn_to_db_refs(urn3)
    assert(db_refs_3 == {'MESH': 'C0012144'})

    # agi-nicmcelltype
    urn4 = 'urn:agi-ncimcelltype:C0242633'
    db_refs_4 = urn_to_db_refs(urn4)
    assert(db_refs_4 == {'MESH': 'C0242633'})

    # agi-meshdist
    urn5 = 'urn:agi-meshdis:Paramyotonia%20Congenita'
    db_refs_5 = urn_to_db_refs(urn5)
    assert(db_refs_5 == {'MESH': 'Paramyotonia%20Congenita'})

    # agi-gocomplex
    urn6 = 'urn:agi-gocomplex:0005610'
    db_refs_6 = urn_to_db_refs(urn6)
    assert(db_refs_6 == {'GO': '0005610'})

    # agi-go
    urn7 = 'urn:agi-go:0001515'
    db_refs_7 = urn_to_db_refs(urn7)
    assert(db_refs_7 == {'GO': '0001515'})

    #agi-ncimtissue
    urn8 = 'urn:agi-ncimtissue:C0007807'
    db_refs_8 = urn_to_db_refs(urn8)
    assert(db_refs_8 == {'MESH': 'C0007807'})
