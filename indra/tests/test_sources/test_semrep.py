import os
from indra.sources import semrep
from indra.statements import *

HERE = os.path.dirname(os.path.abspath(__file__))


def test_semrep_xml():
    fname = os.path.join(HERE, 'resources', 'semrep_example1.xml')
    sp = semrep.process_xml_file(fname)
    assert len(sp.statements) == 3, sp.statements
    assert isinstance(sp.statements[0], Activation)
    assert sp.statements[0].subj.db_refs['UMLS'] == 'C3192263'
    sp = semrep.process_xml_file(fname, use_gilda_grounding=True)
    assert len(sp.statements) == 3, sp.statements
    assert sp.statements[0].subj.db_refs.get('CHEBI') == 'CHEBI:63637'
