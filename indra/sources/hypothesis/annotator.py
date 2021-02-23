__all__ = ['statement_to_annotations']

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

    if evidence.text_refs.get('PMCID'):
        url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/%s/' % \
            evidence.text_refs['PMCID']
    elif evidence.pmid:
        url = 'https://pubmed.ncbi.nlm.nih.gov/%s/' % evidence.pmid
    elif evidence.text_refs.get('URL'):
        url = evidence.text_refs['URL']
    else:
        return None
    return {
        'url': url,
        'target_text': evidence.text,
        'tags': [evidence.source_api]
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
            db_ns, db_id = get_grounding(agent_wc.db_refs,
                                         grounding_ns)
            if not db_ns:
                continue
            identifiers_url = \
                identifiers.get_identifiers_url(db_ns, db_id)
            grounding_text = '[%s](%s)' % (agent_wc.name, identifiers_url)
            insert_len = len(grounding_text) - agent_wc.coords[1] + \
                agent_wc.coords[0]
            inserts.append((agent_wc.coords[0], insert_len))
            before_part = annotation_text[:agent_wc.coords[0]]
            after_part = annotation_text[agent_wc.coords[1]:]
            annotation_text = ''.join([before_part, grounding_text,
                                       after_part])
    return annotation_text
