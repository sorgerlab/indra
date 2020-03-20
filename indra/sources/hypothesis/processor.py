import re
import logging
from indra.statements import BioContext, RefContext
from indra.preassembler.grounding_mapper.standardize import \
    standardize_db_refs, name_from_grounding

logger = logging.getLogger(__name__)


class HypothesisProcessor:
    def __init__(self, annotations, reader=None, grounder=None):
        self.annotations = annotations
        self.statements = []
        if reader is None:
            from indra.sources import reach
            self.reader = reach.process_text
        if grounder is None:
            from gilda import ground
            self.grounder = ground

    def extract_statements(self):
        for annotation in self.annotations:
            stmts = self.stmts_from_annotation(annotation)
            if stmts:
                self.statements += stmts

    def stmts_from_annotation(self, annotation):
        text = annotation.get('text')
        if not text:
            return []
        parts = [t for t in text.split('\n') if t]
        text = parts[0]
        rp = self.reader(text)
        if not rp or not rp.statements:
            logger.warning('Could not extract any statements from %s'
                           % text)
            return []

        contexts = {}
        # We assume that all other parts are related to context
        for part in parts[1:]:
            context_dict = get_context_entry(part, self.grounder, text)
            if context_dict:
                contexts.update(context_dict)
        bio_context = BioContext(**contexts) if contexts else None
        text_refs = get_text_refs(annotation['uri'])
        # In case we got multiple statements out, we apply the same
        # annotations to each
        for stmt in rp.statements:
            # There is expected to be exactly one evidence in all cases
            # but this is still a good way to work with it
            for ev in stmt.evidence:
                ev.source_api = 'hypothes.is'
                ev.text = text
                ev.text_refs = text_refs
                if 'PMID' in text_refs:
                    ev.pmid = text_refs['PMID']
                ev.annotations['hypothes.is'] = annotation
                ev.context = bio_context
        return rp.statements


def get_context_entry(entry, grounder, sentence):
    match = re.match(r'(.*): (.*)', entry)
    if not match:
        return None
    context_type, context_txt = match.groups()
    if context_type not in allowed_contexts:
        logger.warning('Unknown context type %s' % context_type)
        return None

    terms = grounder(context_txt, context=sentence)
    if not terms:
        logger.warning('Could not ground %s context: %s'
                       % (context_type, context_txt))
    db_refs = {}
    if terms:
        db_refs = standardize_db_refs({terms[0].term.db:
                                       terms[0].term.id})
    db_refs['TEXT'] = context_txt
    standard_name = None
    if terms:
        standard_name = name_from_grounding(terms[0].term.db,
                                            terms[0].term.id)
    name = standard_name if standard_name else context_txt
    context = RefContext(name=name, db_refs=db_refs)
    return {context_type: context}


def get_text_refs(url):
    text_refs = {'URL': url}
    match = re.match(r'https://www.ncbi.nlm.nih.gov/pubmed/(\d+)', url)
    if match:
        text_refs['PMID'] = match.groups()[0]
    match = re.match(r'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC(\d+)/',
                     url)
    if match:
        text_refs['PMCID'] = match.groups()[0]
    return text_refs


allowed_contexts = {
    'Location': 'location',
    'Cell line': 'cell_line',
    'Cell type': 'cell_type',
    'Organ': 'organ',
    'Disease': 'disease',
    'Species': 'species'
}
