from indra.assemblers.english import EnglishAssembler
from indra.databases import identifiers
from indra.statements.agent import get_grounding, default_ns_order


grounding_ns = default_ns_order + \
    ['NCIT', 'PUBCHEM', 'CHEMBL']


def statement_to_annotations(stmt, annotate_agents=True):
    annotation_text = get_annotation_text(stmt,
                                          annotate_agents=annotate_agents)
    annotations = []
    for ev in stmt.evidence:
        annot = evidence_to_annotation(ev)
        if annot is None:
            continue
        annot['annotation'] = annotation_text
        annotations.append(annot)
    return annotations


def evidence_to_annotation(evidence):
    if not evidence.text:
        return None

    if 'PMCID' in evidence.text_refs:
        url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/%s/' % \
            evidence.text_refs['PMCID']
    elif evidence.pmid:
        url = 'https://pubmed.ncbi.nlm.nih.gov/%s/' % evidence.pmid
    else:
        return None
    return {
        'url': url,
        'target_text': evidence.text
    }


def get_annotation_text(stmt, annotate_agents=True):
    ea = EnglishAssembler(stmts=[stmt])
    annotation_text = ea.make_model()
    if annotate_agents:
        inserts = []
        for agent_wc in ea.stmt_agents[0]:
            for insert_begin, insert_len in inserts:
                if insert_begin < agent_wc.coords[0]:
                    agent_wc.update_coords(insert_len)
            db_ns, db_id = get_grounding(agent_wc.db_refs)
            if not db_ns:
                continue
            grounding_text = '[%s:%s]' % (db_ns, db_id)
            inserts.append((agent_wc.coords[1], len(grounding_text)))
            before_part = annotation_text[:agent_wc.coords[1]]
            after_part = annotation_text[agent_wc.coords[1]:]
            annotation_text = ''.join([before_part, grounding_text,
                                       after_part])
    return annotation_text
