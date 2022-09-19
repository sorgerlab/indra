import csv
from dataclasses import dataclass, field
from io import TextIOWrapper
import logging
import tarfile
from typing import Any, Dict
import tqdm
from indra.ontology.standardize import get_standard_agent
from indra.statements import *

logger = logging.getLogger(__name__)

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
            # FIXME: look at refined polarity and flip Statement polarity
            # if necessary
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
        return

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
            self.ann_file = TextIOWrapper(fh.extractfile('%s_%s.ann' % key),
                                          encoding='utf-8')
            self.txt_file = TextIOWrapper(fh.extractfile('%s_%s.txt' % key),
                                          encoding='utf-8')
            self.text_lines = self.txt_file.readlines()
            self.elements = process_annotations(self.ann_file)
        self.line_offsets = [0]
        for idx, line in enumerate(self.text_lines[:-1]):
            self.line_offsets.append(self.line_offsets[idx] + len(line))

    def get_sentence_for_offset(self, offset):
        for idx in range(len(self.line_offsets)):
            if self.line_offsets[idx+1] > offset:
                return self.text_lines[idx].strip()

    def find_regulation(self):
        pass


def process_annotations(ann_file):
    elements = {}
    reader = csv.reader(ann_file, delimiter='\t')
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
        elif parts[0] in standoff_event_types:
            event = Event(uid, parts[0], int(parts[1]), int(parts[2]), value)
            elements[uid] = event
        elif parts[0] == 'Confidence':
            if len(row) == 3:
                elements[parts[1]].confidence_val = float(value)
            else:
                elements[parts[1]].confidence_level = parts[2]
        elif len(row) == 2:
            if ':' in parts[0]:
                event_type, parent_id = parts[0].split(':')
                event = elements[parent_id]
                assert event_type == event.event_type
            elif parts[0] in {'Subunit-Complex', 'Protein-Component'}:
                event_type = parts[0]
                event = None
            else:
                assert False, row

            arguments = {}
            for element in parts[1:]:
                role, arg_uid = element.split(':')

                # Some regulations are defined out of order, we need a
                # placeholder for these elements that can be resolved later
                element_obj = elements.get(arg_uid, Unresolved(arg_uid))

                if role in arguments:
                    if not isinstance(arguments[role], list):
                        arguments[role] = [arguments[role]]
                    arguments[role].append(element_obj)
                else:
                    arguments[role] = element_obj
            regulation = Regulation(uid, event, arguments)
            elements[uid] = regulation
        else:
            assert False, row

    # We now need to resolve Unresolved regulation references. At this point
    # it's enough if we take them from the elements dict since they would
    # now be resolved at that level.
    for uid, element in elements.items():
        if isinstance(element, Regulation):
            if isinstance(element.event, Unresolved):
                element.event = elements[element.event.uid]
            for k, v in element.arguments.items():
                if isinstance(v, Unresolved):
                    element.arguments[k] = elements[v.uid]
    return elements


@dataclass
class Entity:
    uid: str
    entity_type: str
    start: int
    end: int
    text: str
    references: Dict[str, str] = field(default_factory=dict)


@dataclass
class Event:
    uid: str
    event_type: str
    start: int
    end: int
    text: str


@dataclass
class Regulation:
    uid: str
    event: Event
    arguments: Dict[str, Any]
    confidence_val: float = None
    confidence_level: str = None


@dataclass
class Unresolved:
    uid: str


standoff_event_types = {
    'Binding',
    'Acetylation',
    'Phosphorylation',
    'Dephosphorylation',
    'Methylation',
    'Ubiquitination',
    'Regulation',
    'Positive_regulation',
    'Negative_regulation',
    'Gene_expression',
    'Transcription',
    'Localization',
}