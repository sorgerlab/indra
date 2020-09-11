import os
import json
from indra.sources import biofactoid

here = os.path.dirname(os.path.abspath(__file__))


def test_process_document():
    doc_json = os.path.join(here, 'biofactoid_doc.json')
    with open(doc_json, 'r') as fh:
        doc = json.load(fh)
    bp = biofactoid.process_json([doc])
    assert len(bp.statements) == 2
    assert {s.__class__.__name__ for s in bp.statements} == \
        {'Inhibition', 'Phosphorylation'}
    for stmt in bp.statements:
        agents = stmt.agent_list()
        assert agents[0].name == 'AKT1'
        assert agents[0].db_refs == \
            {'EGID': '207', 'HGNC': '391', 'ENSEMBL': 'ENSG00000142208'}
        assert agents[1].name == 'FOXO3', agents
        assert agents[1].db_refs == \
            {'EGID': '2309', 'HGNC': '3821', 'ENSEMBL': 'ENSG00000118689'}
        ev = stmt.evidence[0]
        assert ev.pmid == '29886111'
        assert ev.text == 'AKT1 inhibits FOXO3 via phosphorylation.'
        assert ev.annotations == \
            {"biofactoid_document": "3d2a77ba-55c1-463b-aff3-9acaa0307b62",
             "created_date": "2020-09-07T16:38:08.837Z",
             "lsatEditedDate": "2020-09-10T01:13:41.528Z"}
        assert ev.text_refs == \
               {'PMID': '29886111',
                'DOI': '10.1016/j.cels.2018.05.004',
                'PII': 'S2405-4712(18)30192-3',
                'PMCID': 'PMC6322215'}
        assert ev.source_api == 'biofactoid'