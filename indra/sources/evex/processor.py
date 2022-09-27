import csv
from collections import defaultdict
from dataclasses import dataclass, field
from io import TextIOWrapper
import itertools
import logging
import tarfile
from typing import Any, Dict
import networkx
import tqdm
from indra.ontology.standardize import get_standard_agent
from indra.statements import *

logger = logging.getLogger(__name__)


class EvexProcessor:
    def __init__(self, relations_table, articles_table, standoff_index):
        self.relations_table = relations_table
        self.articles_table = articles_table
        self.article_lookup = self.build_article_lookup()
        self.standoff_index = standoff_index
        self.statements = []
        self.standoff_cache = {}
        self.provenance_stats = {
            'file_missing': set(),
            'reg_number_not_one': {},
        }

    def print_provenance_stats(self):
        from collections import Counter
        print('Standoff files not found for %s relations' %
              len(self.provenance_stats['file_missing']))
        print('Regulation numbers other than one found:')
        reg_nums = Counter(
            self.provenance_stats['reg_number_not_one'].values()).most_common()
        for num, cnt in reg_nums:
            print('- %s: %s' % (num, cnt))

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table),
                             desc='Processing Evex relations'):
            pol_idx = 1 if row.refined_polarity == 'Negative' else 0
            stmt_types = type_indra_mappings.get(row.refined_type)
            if not stmt_types:
                continue
            stmt_type = stmt_types[pol_idx]
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

            standoff = self.get_standoff_for_event(row, article_prefix,
                                                   article_id)
            if not standoff:
                continue

            standoff_event_types = type_standoff_mappings[row.refined_type]

            regs = standoff.find_potential_regulations(source_id, target_id)
            matching_regs = set()
            for reg in regs:
                source_paths = reg.paths_to_entrez_id(source_id)
                source_annotated_paths = [standoff.annotate_path(source_path)
                                          for source_path in source_paths]
                target_paths = reg.paths_to_entrez_id(target_id)
                target_annotated_paths = [standoff.annotate_path(target_path)
                                          for target_path in target_paths]
                if row.refined_polarity in {'Positive', 'Unspecified'}:
                    constraints = [
                        {
                            'source': {'Positive_regulation', 'Regulation'},
                            'target': {standoff_event_types[0]}
                         },
                        {
                            'source': {'Negative_regulation'},
                            'target': {standoff_event_types[1]}
                        }
                    ]
                else:
                    constraints = [
                        {
                            'source': {'Positive_regulation', 'Regulation'},
                            'target': {standoff_event_types[1]}
                        },
                        {
                            'source': {'Negative_regulation'},
                            'target': {standoff_event_types[0]}
                        }
                    ]
                for source_path, target_path in \
                        itertools.product(source_annotated_paths,
                                          target_annotated_paths):
                    for constraint in constraints:
                        if source_path[0] in constraint['source'] \
                                and source_path[1] == 'Cause' \
                                and target_path[1] == 'Theme' \
                                and constraint['target'] in target_path:
                            matching_regs.add(reg)

            if not matching_regs:
                label = '\n'.join(['%s: %s' % (k, getattr(row, k))
                                   for k in self.relations_table.columns])
                standoff.save_potential_regulations(source_id, target_id,
                                                    key=standoff.key,
                                                    label=label)
                texts = []
            else:
                texts = [standoff.get_sentence_for_offset(reg.event.start)
                         for reg in matching_regs]

            for text in texts:
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
            self.provenance_stats['file_missing'].add(key)
            return None
        standoff = EvexStandoff(standoff_file, key)
        self.standoff_cache[key] = standoff
        return standoff

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
        self.key = key
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

    def find_exact_regulations(self, cause_entrez_id, theme_entrez_id):
        regs = []
        for uid, element in self.elements.items():
            if isinstance(element, Regulation):
                if {cause_entrez_id, theme_entrez_id} == element.entrez_ids:
                    regs.append(element)
        return regs

    def find_potential_regulations(self, cause_entrez_id, theme_entrez_id):
        regs = []
        for uid, element in self.elements.items():
            if isinstance(element, Regulation):
                if {cause_entrez_id, theme_entrez_id} <= element.entrez_ids:
                    regs.append(element)
        return regs

    def save_potential_regulations(self, cause_entrez_id, theme_entrez_id, key,
                                   label):
        import pystow
        file_key = '_'.join(list(key) + [cause_entrez_id, theme_entrez_id])
        fname = pystow.join('evex', 'debug', name='%s.pdf' % file_key)
        regs = self.find_potential_regulations(cause_entrez_id, theme_entrez_id)
        if not regs:
            return []
        graph = networkx.compose_all([reg.graph for reg in regs])
        ag = networkx.nx_agraph.to_agraph(graph)
        ag.graph_attr['label'] = label
        ag.draw(fname, prog='dot')
        return regs

    def annotate_path(self, path_nodes):
        root = self.elements[path_nodes[0]]
        path_info = [root.event.get_type()]
        for source, target in zip(path_nodes[:-1], path_nodes[1:]):
            path_info.append(root.graph.edges[(source, target)]['label'])
            if isinstance(self.elements[target], Regulation):
                path_info.append(self.elements[target].event.event_type)
        return path_info


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
            # These events don't have actual objects associated with them so
            # we create placeholder events just to propagate the type
            elif parts[0] in {'Subunit-Complex', 'Protein-Component'}:
                event_type = parts[0]
                event = PlaceholderEvent(event_type)
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

    # Now that everything is resolved, we can initialize the regulations
    for uid, element in elements.items():
        if isinstance(element, Regulation):
            element.initialize()

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
class PlaceholderEvent:
    event_type: str

    def get_type(self):
        return self.event_type


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
    # Dynamically created attributes
    entrez_ids = None
    entrez_uid_mappings = None
    graph = None

    def initialize(self):
        # Note this can't be simply post init because of unresolved child
        # objects upon initialization
        self.entrez_ids = self.find_entrez_ids()
        self.entrez_uid_mappings = self.get_entrez_uid_mappings()
        self.graph = self.to_graph()

    def to_graph(self):
        g = networkx.DiGraph()
        add_subgraph(g, self)
        return g

    def draw(self, fname):
        ag = networkx.nx_agraph.to_agraph(self.graph)
        ag.draw(fname, prog='dot')

    def find_entrez_ids(self):
        """Return all Entrez IDs under this regulation."""
        entrez_ids = set()
        for k, v in self.arguments.items():
            if isinstance(v, Regulation):
                entrez_ids |= v.find_entrez_ids()
            elif isinstance(v, Entity):
                entrez_id = v.references.get('EG')
                if entrez_id:
                    entrez_ids.add(entrez_id)
        return entrez_ids

    def get_entrez_uid_mappings(self):
        """Return mappings from Entrez IDs to the UIDs of the nodes where it
        appears."""
        uid_mappings = defaultdict(list)
        for arg_type, arg in self.arguments.items():
            for single_arg in (arg if isinstance(arg, list) else [arg]):
                if isinstance(single_arg, Regulation):
                    for child_entrez_id, child_uids in \
                            single_arg.get_entrez_uid_mappings().items():
                        for child_uid in child_uids:
                            uid_mappings[child_entrez_id].append(child_uid)
                elif isinstance(single_arg, Entity):
                    entrez_id = single_arg.references.get('EG')
                    if entrez_id:
                        uid_mappings[entrez_id].append(single_arg.uid)
        return dict(uid_mappings)

    def paths_to_entrez_id(self, entrez_id):
        """Find a path from the root to a given Entrez ID."""
        uids = self.entrez_uid_mappings.get(entrez_id)
        paths = []
        for uid in uids:
            path_nodes = networkx.shortest_path(self.graph, self.uid, uid)
            if len(path_nodes) == 1:
                breakpoint()
            paths.append(path_nodes)
        return paths


def add_subgraph(g, obj):
    label = '{ID | %s} | {event_type | %s}' % (obj.uid, obj.event.get_type())
    if obj.negation:
        label += '| {negated | %s}' % True
    g.add_node(obj.uid, type='Regulation',
               shape='record',
               label=label)
    for k, v in obj.arguments.items():
        for vv in (v if isinstance(v, list) else [v]):
            if isinstance(vv, Regulation):
                add_subgraph(g, vv)
            else:
                label = '{ID | %s} | {type | %s} | {text | %s}' % \
                    (vv.uid, vv.get_type(), vv.text)
                if isinstance(vv, Entity):
                    egid = vv.references.get('EG')
                    if egid:
                        label += '| {entrez_id | %s}' % egid
                g.add_node(vv.uid,
                           shape='record',
                           label=label)
            g.add_edge(obj.uid, vv.uid, label=k)


@dataclass
class Unresolved:
    uid: str


# The set of event types used in the standoff format
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

# Mapping network relation types to regulation types used in the standoff files
# as well as the one with opposite polarity.
type_standoff_mappings = {
    'Binding': ('Binding', 'Binding'),
    'Catalysis of DNA methylation': ('Methylation', 'Demethylation'),
    'Catalysis of acetylation': ('Acetylation', 'Deacethylation'),
    'Catalysis of glycosylation': ('Glycosylation', 'Deglycosylation'),
    'Catalysis of hydroxylation': ('Hydroxylation', 'Dehydroxylation'),
    'Catalysis of methylation': ('DNA_methylation', 'DNA_demethylation'),
    'Catalysis of phosphorylation': ('Phosphorylation', 'Dephosphorylation'),
    'Catalysis of ubiquitination': ('Ubiquitination', 'Deubiquitination'),
    'Indirect_catalysis of acetylation': ('Acetylation', 'Deacethylation'),
    'Indirect_catalysis of methylation': ('Methylation', 'Demethylation'),
    'Indirect_catalysis of ubiquitination': ('Ubiquitination',
                                             'Deubiquitination'),
    'Indirect_regulation': (None, None),
    'Indirect_regulation of binding': ('Binding', 'Binding'),
    'Indirect_regulation of catabolism': ('Protein_catabolism',
                                          'Protein_catabolism'),
    'Indirect_regulation of expression': ('Gene_expression', 'Gene_expression'),
    'Indirect_regulation of localization': ('Localization', 'Localization'),
    'Indirect_regulation of phosphorylation': ('Phosphorylation',
                                               'Dephosphorylation'),
    'Indirect_regulation of transcription': ('Transcription', 'Transcription'),
    'Regulation': (None, None),
    'Regulation of binding': ('Binding', 'Binding'),
    'Regulation of catabolism': ('Protein_catabolism', 'Protein_catabolism'),
    'Regulation of expression': ('Gene_expression', 'Gene_expression'),
    'Regulation of localization': ('Localization', 'Localization'),
    'Regulation of phosphorylation': ('Phosphorylation', 'Dephosphorylation'),
    'Regulation of transcription': ('Transcription', 'Transcription'),
}

# Network relation type mappings to INDRA Statement types
type_indra_mappings = {
    'Binding': (Complex, Complex),
    'Catalysis of DNA methylation': (Methylation, Demethylation),
    'Catalysis of acetylation': (Acetylation, Deacetylation),
    'Catalysis of glycosylation': (Glycosylation, Deglycosylation),
    'Catalysis of hydroxylation': (Hydroxylation, Dehydroxylation),
    'Catalysis of methylation': (Methylation, Demethylation),
    'Catalysis of phosphorylation': (Phosphorylation, Dephosphorylation),
    'Catalysis of ubiquitination': (Ubiquitination, Deubiquitination),
    'Indirect_catalysis of acetylation': (Acetylation, Deacetylation),
    'Indirect_catalysis of methylation': (Methylation, Demethylation),
    'Indirect_catalysis of ubiquitination': (Ubiquitination, Deubiquitination),
    'Indirect_regulation': None,
    'Indirect_regulation of binding': None,
    'Indirect_regulation of catabolism': None,
    'Indirect_regulation of expression': (IncreaseAmount, DecreaseAmount),
    'Indirect_regulation of localization': None,
    'Indirect_regulation of phosphorylation': (Phosphorylation,
                                               Dephosphorylation),
    'Indirect_regulation of transcription': (IncreaseAmount, DecreaseAmount),
    'Regulation': None,
    'Regulation of binding': None,
    'Regulation of catabolism': None,
    'Regulation of expression': (IncreaseAmount, DecreaseAmount),
    'Regulation of localization': None,
    'Regulation of phosphorylation': (Phosphorylation, Dephosphorylation),
    'Regulation of transcription': (IncreaseAmount, DecreaseAmount),
}
