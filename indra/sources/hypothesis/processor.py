import re
from indra.statements import BioContext, RefContext
from indra.preassembler.grounding_mapper.standardize import \
    standardize_db_refs, name_from_grounding


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
            if stmt:
                self.statements.append(stmt)

    def stmts_from_annotation(self, annotation):
        text = annotation.get('text')
        if not text:
            return []
        parts = [t for t in text.split('\n') if t]
        rp = self.reader(parts[0])
        for part in parts[1:]:
            match = re.match(r'Context: (.*)', part)
            if not match:
                continue
            context_txt = match.groups()[0]
            terms = self.grounder(context_txt)
            db_refs = standardize_db_refs({terms[0].db: terms[0].id}) \
                if terms else {}
            db_refs['TEXT'] = context_txt
            standard_name = name_from_grounding(terms[0].db, terms[0].id)
            name = standard_name if standard_name else context_txt
            # TODO: how can we tell what kind of BioContext this is exactly?
            # Disease, organ, cell type, etc?
            context = RefContext(name=name, db_refs=db_refs)
        for stmt in rp.statements:
            stmt.evidence[0].annotations.extend(annotation)

