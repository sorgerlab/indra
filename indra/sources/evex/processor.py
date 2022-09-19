import tarfile
import tqdm
import csv
from indra.statements import *
from indra.ontology.standardize import get_standard_agent

type_mappings = {
    'Binding': Complex,
    'Catalysis of DNA methylation': Methylation,
    'Catalysis of acetylation': Acetylation,
    'Catalysis of glycosylation': Glycosylation,
    'Catalysis of hydroxylation': Hydroxylation,
    'Catalysis of methylation': Methylation,
    'Catalysis of phosphorylation': Phosphorylation,
    'Catalysis of ubiquitination': Ubiquitination,
    'Indirect_catalysis of acetylation': Acetylation,
    'Indirect_catalysis of methylation': Methylation,
    'Indirect_catalysis of ubiquitination': Ubiquitination,
    'Indirect_regulation': None,
    'Indirect_regulation of binding': None,
    'Indirect_regulation of catabolism': None,
    'Indirect_regulation of expression': IncreaseAmount,
    'Indirect_regulation of localization': None,
    'Indirect_regulation of phosphorylation': Phosphorylation,
    'Indirect_regulation of transcription': IncreaseAmount,
    'Regulation': None,
    'Regulation of binding': None,
    'Regulation of catabolism': None,
    'Regulation of expression': IncreaseAmount,
    'Regulation of localization': None,
    'Regulation of phosphorylation': Phosphorylation,
    'Regulation of transcription': IncreaseAmount,
}


class EvexProcessor:
    def __init__(self, relations_table, articles_table, standoff_index):
        self.relations_table = relations_table
        self.articles_table = articles_table
        self.article_lookup = self.build_article_lookup()
        self.standoff_index = standoff_index
        self.statements = []

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table),
                             desc='Processing Evex relations'):
            stmt_type = type_mappings.get(row.refined_type)
            if not stmt_type:
                continue
            subj_agent = get_standard_agent(
                'EGID:%s' % str(row.source_entrezgene_id),
                db_refs={'EGID': str(row.source_entrezgene_id)}
            )
            obj_agent = get_standard_agent(
                'EGID:%s' % str(row.target_entrezgene_id),
                db_refs={'EGID': str(row.target_entrezgene_id)}
            )

            article_prefix, article_id = \
                self.article_lookup.get(row.general_event_id)

            text_refs = {article_prefix: article_id}
            pmid = article_id if article_prefix == 'PMID' else None

            text = self.get_text_for_event(row.general_event_id,
                                           article_prefix,
                                           article_id)

            ev = Evidence(source_api='evex',
                          pmid=pmid,
                          text_refs=text_refs)

            if stmt_type == Complex:
                stmt = Complex([subj_agent, obj_agent], evidence=[ev])
            else:
                stmt = stmt_type(subj_agent, obj_agent, evidence=[ev])
            self.statements.append(stmt)

    def get_text_for_event(self, event_id, article_prefix, article_id):
        key = (
            'pmc' if article_prefix == 'PMCID' else 'pubmed',
            article_id[3:] if article_prefix == 'PMCID' else article_id
        )
        standoff_file = self.standoff_index.get(key)
        if not standoff_file:
            return None
        else:
            breakpoint()
        standoff = EvexStandoff(standoff_file, key)


    def build_article_lookup(self):
        article_lookup = {}
        for row in self.articles_table.itertuples():
            prefix, article_id = row.article_id.split(': ')
            if prefix == 'PMCID':
                if not article_id.startswith('PMC'):
                    article_id = 'PMC' + article_id
                article_lookup[row.general_event_id] = ('PMCID', article_id)
            elif prefix == 'PMID':
                article_lookup[row.general_event_id] = ('PMID', article_id)
            else:
                ValueError('Unexpected article type: %s' % prefix)
        return article_lookup


class EvexStandoff:
    def __init__(self, standoff_file, key):
        with tarfile.open(standoff_file, 'r:gz') as fh:
            self.ann_file = fh.extractfile('%s_%s.ann' % key)
            self.txt_file = fh.extractfile('%s_%s.txt' % key)
            self.elements = process_annotations(self.ann_file)


def process_annotations(ann_file):
    elements = {}
    import io
    decoded = io.TextIOWrapper(ann_file, encoding="utf-8")
    reader = csv.reader(decoded, delimiter='\t')
    for row in reader:
        uid = row[0]
        assert len(row) == 2 or len(row) == 3
        value = row[2] if len(row) == 3 else None
        parts = row[1].split()
        if parts[0] in {'GGP', 'Entity'}:
            entity = Entity(uid, parts[0], int(parts[1]), int(parts[2]), value)
            elements[uid] = entity
        elif parts[0] == 'Reference':
            ref_ns, ref_id = parts[2].split(':', maxsplit=1)
            elements[parts[1]].references[ref_ns] = ref_id
        elif parts[0] in standoff_relation_types:
            relation = Relation(uid, parts[0], int(parts[1]), int(parts[2]),
                                value)
            elements[uid] = relation
        elif parts[0] == 'Confidence':
            if len(row) == 3:
                elements[parts[1]].confidence_val = float(value)
            else:
                elements[parts[1]].confidence_level = parts[1]
        elif len(row) == 2:
            if ':' in parts[0]:
                relation_type, parent_id = parts[0].split(':')
                parent = elements[parent_id]
            elif parts[0] in {'Subunit-Complex', 'Protein-Component'}:
                relation_type = parts[0]
                parent = None
            else:
                assert False, row

            arguments = {}
            for element in parts[1:]:
                role, arg_uid = element.split(':')
                if role in arguments:
                    if not isinstance(arguments[role], list):
                        arguments[role] = [arguments[role][0]]
                    arguments[role].append(elements[arg_uid])
            regulation = Regulation(uid, parent, arguments)
            elements[uid] = regulation
        else:
            assert False, row
    return elements


class Entity:
    def __init__(self, uid, entity_type, start, end, text, references=None):
        self.uid = uid
        self.type = entity_type
        self.start = start
        self.end = end
        self.text = text
        self.references = references if references else {}


class Relation:
    def __init__(self, uid, relation_type, start, end, text):
        self.uid = uid
        self.relation_type = relation_type
        self.start = start
        self.end = end
        self.text = text


class Regulation:
    def __init__(self, uid, parent, arguments, confidence_val=None,
                 confidence_level=None):
        self.uid = uid
        self.parent = parent
        self.arguments = arguments
        self.confidence_val = confidence_val
        self.confidence_level = confidence_level


standoff_relation_types = {
    'Binding',
    'Phosphorylation',
    'Dephosphorylation',
    'Methylation',
    'Regulation',
    'Positive_regulation',
    'Negative_regulation',
    'Gene_expression',
    'Transcription',
    'Localization',
}