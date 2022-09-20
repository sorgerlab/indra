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
        self.standoff_cache = {}

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table),
                             desc='Processing Evex relations'):
            # FIXME: look at refined polarity and flip Statement polarity
            # if necessary
            stmt_type = type_mappings.get(row.refined_type)
            if not stmt_type:
                continue
            source_id = str(row.source_entrezgene_id)
            target_id = str(row.target_entrezgene_id)
            subj_agent = get_standard_agent('EGID:%s' % source_id,
                                            db_refs={'EGID': source_id})
            obj_agent = get_standard_agent('EGID:%s' % target_id,
                                           db_refs={'EGID': target_id})

            article_prefix, article_id = \
                self.article_lookup.get(row.general_event_id)

            text_refs = {article_prefix: article_id}
            pmid = article_id if article_prefix == 'PMID' else None

            st = self.get_standoff_for_event(row, article_prefix, article_id)
            if st:
                text = self.get_text_from_standoff(st, source_id, target_id)
            else:
                text = None

            ev = Evidence(source_api='evex',
                          pmid=pmid,
                          text_refs=text_refs,
                          text=text)

            if stmt_type == Complex:
                stmt = Complex([subj_agent, obj_agent], evidence=[ev])
            else:
                stmt = stmt_type(subj_agent, obj_agent, evidence=[ev])
            self.statements.append(stmt)

    def get_standoff_for_event(self, row, article_prefix, article_id):
        key = (
            'pmc' if article_prefix == 'PMCID' else 'pubmed',
            article_id[3:] if article_prefix == 'PMCID' else article_id
        )
        if key in self.standoff_cache:
            return self.standoff_cache[key]
        standoff_file = self.standoff_index.get(key)
        if not standoff_file:
            return None
        standoff = EvexStandoff(standoff_file, key)
        self.standoff_cache[key] = standoff
        return standoff

    @staticmethod
    def get_text_from_standoff(standoff, source_id, target_id):
        regs = standoff.find_regulations(source_id, target_id)
        # FIXME: implement more sophisticated matching
        if len(regs) != 1:
            return None
        return standoff.get_sentence_for_offset(regs[0].event.start)

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
        for idx in range(len(self.line_offsets) - 1):
            if self.line_offsets[idx+1] > offset:
                return self.text_lines[idx].strip()
        return len(self.line_offsets) - 1

    def find_regulations(self, cause_entrez_id, theme_entrez_id):
        regs = []
        for uid, element in self.elements.items():
            if isinstance(element, Regulation):
                entrez_ids = element.find_entrez_ids()
                if {cause_entrez_id, theme_entrez_id} == entrez_ids:
                    regs.append(element)
        return regs


def process_annotations(ann_file):
    elements = {}
    reader = csv.reader(ann_file, delimiter='\t', quotechar=None)
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
            # Negation confidence
            if isinstance(parts[1], Negation):
                elements[parts[1]].confidence = float(value)
            # Regulation confidence value
            elif len(row) == 3:
                elements[parts[1]].confidence_val = float(value)
            # Regulation confidence level
            else:
                elements[parts[1]].confidence_level = parts[2]
        elif parts[0] == 'Negation':
            elements[uid] = Negation(uid)
            elements[parts[1]].negation = elements[uid]
        elif parts[0] == 'Speculation':
            elements[uid] = Speculation(uid)
            elements[parts[1]].speculation = elements[uid]
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
            print(row)
            break

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
class Negation:
    uid: str
    confidence: float = None


@dataclass
class Speculation:
    uid: str
    confidence: float = None


@dataclass
class Entity:
    uid: str
    entity_type: str
    start: int
    end: int
    text: str
    references: Dict[str, str] = field(default_factory=dict)

    def get_type(self):
        return self.entity_type


@dataclass
class Event:
    uid: str
    event_type: str
    start: int
    end: int
    text: str

    def get_type(self):
        return self.event_type


@dataclass
class Regulation:
    uid: str
    event: Event
    arguments: Dict[str, Any]
    confidence_val: float = None
    confidence_level: str = None
    negation: Negation = None
    speculation: Speculation = None

    def to_graph(self):
        from networkx import DiGraph
        g = DiGraph()
        add_subgraph(g, self)
        return g

    def draw(self, fname):
        import networkx
        ag = networkx.nx_agraph.to_agraph(self.to_graph())
        ag.draw(fname, prog='dot')

    def find_entrez_ids(self):
        entrez_ids = set()
        for k, v in self.arguments.items():
            if isinstance(v, Regulation):
                entrez_ids |= v.find_entrez_ids()
            elif isinstance(v, Entity):
                entrez_id = v.references.get('EG')
                if entrez_id:
                    entrez_ids.add(entrez_id)
        return entrez_ids


def add_subgraph(g, obj):
    g.add_node(obj.uid, type='Regulation')
    for k, v in obj.arguments.items():
        if isinstance(v, Regulation):
            add_subgraph(g, v)
        else:
            g.add_node(v.uid, type=v.get_type())
        g.add_edge(obj.uid, v.uid, type=k)


@dataclass
class Unresolved:
    uid: str


standoff_event_types = {
    'Binding',
    'Acetylation',
    'Deacetylation',
    'Phosphorylation',
    'Dephosphorylation',
    'DNA_methylation',
    'DNA_demethylation',
    'Glycosylation',
    'Deglycosylation',
    'Hydroxylation',
    'Dehydroxylation',
    'Methylation',
    'Demethylation',
    'Ubiquitination',
    'Deubiquitination',
    'Regulation',
    'Positive_regulation',
    'Negative_regulation',
    'Gene_expression',
    'Catalysis',
    'Transcription',
    'Localization',
    'Protein_catabolism',
}