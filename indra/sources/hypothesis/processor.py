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
        rp = self.reader(parts[0])
        if not rp or not rp.statements:
            logger.warning('Could not extract ant statements from %s'
                           % parts[0])
            return []
        contexts = {}
        for part in parts[1:]:
            match = re.match(r'(.*): (.*)', part)
            if not match:
                continue
            context_type, context_txt = match.groups()
            if context_type not in allowed_contexts:
                logger.warning('Unknown context type %s' % context_type)
                continue

            terms = self.grounder(context_txt)
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
            contexts[allowed_contexts[context_type]] = context
        bio_context = BioContext(**contexts) if contexts else None
        for stmt in rp.statements:
            stmt.evidence.source_api = 'hypothes.is'
            stmt.evidence[0].annotations.update(annotation)
            stmt.evidence[0].context = bio_context
        return rp.statements


allowed_contexts = {
    'Location': 'location',
    'Cell line': 'cell_line',
    'Cell type': 'cell_type',
    'Organ': 'organ',
    'Disease': 'disease',
    'Species': 'species'
}
