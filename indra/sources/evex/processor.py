import copy
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
from indra.statements.validate import validate_statement
from indra.statements import *

logger = logging.getLogger(__name__)


class EvexProcessor:
    """A processor to extract INDRA Statements from EVEX relations."""
    def __init__(self, relations_table, articles_table, standoff_index):
        self.relations_table = relations_table
        self.articles_table = articles_table
        # Build an index of
        self.article_lookup = self.build_article_lookup()
        self.standoff_index = standoff_index
        self.statements = []
        self.standoff_cache = {}

    def process_statements(self):
        for row in tqdm.tqdm(self.relations_table.itertuples(),
                             total=len(self.relations_table),
                             desc='Processing Evex relations'):
            self.statements += self.process_row(row)

    def process_row(self, row):
        """Process a row in the relations table into INDRA Statements."""
        pol_idx = 1 if row.refined_polarity == 'Negative' else 0
        stmt_types = type_indra_mappings.get(row.refined_type)
        if not stmt_types:
            return []
        stmt_type = stmt_types[pol_idx]
        source_id = str(row.source_entrezgene_id)
        target_id = str(row.target_entrezgene_id)
        subj_agent = get_standard_agent('EGID:%s' % source_id,
                                        db_refs={'EGID': source_id})
        obj_agent = get_standard_agent('EGID:%s' % target_id,
                                       db_refs={'EGID': target_id})

        article_keys = self.article_lookup.get(row.general_event_id)
        stmts = []
        for article_prefix, article_id in article_keys:
            text_refs = {article_prefix: article_id}
            pmid = article_id if article_prefix == 'PMID' else None

            standoff = self.get_standoff_for_event(article_prefix, article_id)
            if not standoff:
                evidence_info = [{}]
            else:
                evidence_info = find_evidence_info(standoff, source_id,
                                                   target_id, row.refined_type,
                                                   row.refined_polarity)
            for ev_info in evidence_info:
                annotations = {
                    'evex_relation_type': row.refined_type,
                    'evex_polarity': row.refined_polarity,
                    'evex_general_event_id': str(row.general_event_id),
                    'evex_standoff_regulation_id': ev_info.get('regulation_uid')
                }
                if ev_info.get('subj_coords'):
                    annotations['agents'] = \
                        {'coords': [ev_info['subj_coords'],
                                    ev_info['obj_coords']]}
                ev = Evidence(source_api='evex',
                              pmid=pmid,
                              text_refs=text_refs,
                              text=ev_info.get('text'),
                              annotations=annotations)
                subj = copy.deepcopy(subj_agent)
                obj = copy.deepcopy(obj_agent)
                if ev_info.get('subj_text'):
                    subj.db_refs['TEXT'] = ev_info.get('subj_text')
                if ev_info.get('obj_text'):
                    obj.db_refs['TEXT'] = ev_info.get('obj_text')
                if stmt_type == Complex:
                    stmt = Complex([subj, obj], evidence=[ev])
                else:
                    stmt = stmt_type(subj, obj, evidence=[ev])
                validate_statement(stmt)
                stmts.append(stmt)
        return stmts

    def get_standoff_for_event(self, article_prefix, article_id):
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

    def build_article_lookup(self):
        """Build a lookup for articles corresponding to event IDs."""
        article_lookup = defaultdict(list)
        for row in self.articles_table.itertuples():
            prefix, article_id = row.article_id.split(': ')
            if prefix == 'PMCID':
                if not article_id.startswith('PMC'):
                    article_id = 'PMC' + article_id
                article_lookup[row.general_event_id].append(
                    ('PMCID', article_id))
            elif prefix == 'PMID':
                article_lookup[row.general_event_id].append(
                    ('PMID', article_id))
            else:
                ValueError('Unexpected article type: %s' % prefix)
        return dict(article_lookup)


def find_evidence_info(standoff, source_id, target_id, event_type,
                       polarity):
    """Given a standoff, find all regulations matching a relation row
    and return corresponding evidence info."""
    potential_regs = standoff.find_potential_regulations(source_id,
                                                         target_id)
    matching_reg_info = []
    for reg in potential_regs:
        source_paths = reg.paths_to_entrez_id(source_id)
        source_annotated_paths = [standoff.annotate_path(source_path)
                                  for source_path in source_paths]
        target_paths = reg.paths_to_entrez_id(target_id)
        target_annotated_paths = [standoff.annotate_path(target_path)
                                  for target_path in target_paths]

        if event_type == 'Binding':
            for source_path, target_path in \
                    itertools.product(source_annotated_paths,
                                      target_annotated_paths):
                if 'Binding' in source_path and 'Binding' in target_path:
                    source_reg_idx = source_path.index('Binding')
                    target_reg_idx = target_path.index('Binding')
                    if source_path[source_reg_idx + 1] == 'Theme' and \
                            target_path[target_reg_idx + 1] == 'Theme':
                        matching_reg_info.append(
                            get_regulation_info(standoff, reg,
                                                source_path[-1],
                                                target_path[-1]))
        else:
            pos_event_type, neg_event_type = \
                type_standoff_mappings[event_type]
            polarity_is_positive = (polarity in {'Positive', 'Unspecified'})
            constraints = [
                {
                    'source': {'Positive_regulation', 'Regulation',
                               'Catalysis'},
                    'target': (pos_event_type if polarity_is_positive
                               else neg_event_type)
                },
                {
                    'source': {'Negative_regulation'},
                    'target': (neg_event_type if polarity_is_positive
                               else pos_event_type)
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
                        matching_reg_info.append(
                            get_regulation_info(standoff, reg,
                                                source_path[-1],
                                                target_path[-1]))

    if not matching_reg_info:
        if len(potential_regs) == 1:
            txt = standoff.get_sentence_for_offset(
                potential_regs[0].event.start)
            matching_reg_info = [{'text': txt,
                                  'regulation_uid': potential_regs[0].uid}]
        else:
            data = {'source_id': source_id,
                    'target_id': target_id,
                    'event_type': event_type,
                    'polarity': polarity}
            label = '\n'.join(['%s: %s' % (k, v) for k, v in data.items()])
            standoff.save_potential_regulations(source_id, target_id,
                                                key=standoff.key,
                                                label=label)
            matching_reg_info = [{}]

    return matching_reg_info


def get_regulation_info(standoff, regulation, source_uid, target_uid):
    text = standoff.get_sentence_for_offset(regulation.event.start)
    subj = standoff.elements[source_uid]
    subj_text = subj.text
    subj_coord = [standoff.get_sentence_relative_offset(subj.start),
                  standoff.get_sentence_relative_offset(subj.end)]
    obj = standoff.elements[target_uid]
    obj_text = obj.text
    obj_coord = [standoff.get_sentence_relative_offset(obj.start),
                 standoff.get_sentence_relative_offset(obj.end)]
    return {'text': text,
            'subj_text': subj_text,
            'subj_coord': subj_coord,
            'obj_text': obj_text,
            'obj_coord': obj_coord,
            'regulation_uid': regulation.uid}


def get_sentence_for_offset(text_lines, line_offsets, offset):
    """Return a text line for a given offset based on line offsets."""
    for idx in range(len(line_offsets) - 1):
        if line_offsets[idx + 1] > offset:
            return text_lines[idx].strip()
    return text_lines[-1]


def get_sentence_relative_offset(line_offsets, offset):
    """Return an offset relative to the sentence it is in."""
    for idx in range(len(line_offsets) - 1):
        if line_offsets[idx + 1] > offset:
            return offset - line_offsets[idx]
    return offset - line_offsets[-1]


class EvexStandoff:
    """Represent an EVEX standoff file's contents as a set of objects."""
    def __init__(self, standoff_file, key):
        self.key = key
        # We need to get the content of the text lines corresponding to
        # the standoff annotations, and then process the annotations from
        # the annotation file.
        with tarfile.open(standoff_file, 'r:gz') as fh:
            ann_file = TextIOWrapper(fh.extractfile('%s_%s.ann' % key),
                                     encoding='utf-8')
            txt_file = TextIOWrapper(fh.extractfile('%s_%s.txt' % key),
                                     encoding='utf-8')
            self.text_lines = txt_file.readlines()
            self.elements = process_annotations(ann_file)
        # To be able to linearly index into sentences broken up into separate
        # lines, we build an index of line offsets
        self.line_offsets = [0]
        for idx, line in enumerate(self.text_lines[:-1]):
            self.line_offsets.append(self.line_offsets[idx] + len(line))

    def get_sentence_for_offset(self, offset):
        """Return the sentence for a given offset in the standoff annotation."""
        return get_sentence_for_offset(self.text_lines, self.line_offsets,
                                       offset)

    def get_sentence_relative_offset(self, offset):
        """Return an offset relative to the sentence it is in."""
        return get_sentence_relative_offset(self.line_offsets, offset)

    def find_exact_regulations(self, cause_entrez_id, theme_entrez_id):
        """Find regulations that only contain the given entrez IDs."""
        regs = []
        for uid, element in self.elements.items():
            if isinstance(element, Regulation):
                if {cause_entrez_id, theme_entrez_id} == element.entrez_ids:
                    regs.append(element)
        return regs

    def find_potential_regulations(self, cause_entrez_id, theme_entrez_id):
        """Find regulations that contain the given entrez IDs."""
        regs = []
        for uid, element in self.elements.items():
            if isinstance(element, Regulation):
                if {cause_entrez_id, theme_entrez_id} <= element.entrez_ids:
                    regs.append(element)
        return regs

    def save_potential_regulations(self, cause_entrez_id, theme_entrez_id, key,
                                   label):
        """Save potential regulation graphs for review/debugging."""
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
        """Given a raw path of node IDs, create an annotated path.

        The annotated path contains event types for regulation nodes,
        relation types, and IDs of leaf entity nodes.
        """
        root = self.elements[path_nodes[0]]
        path_info = [root.event.get_type()]
        for source, target in zip(path_nodes[:-1], path_nodes[1:]):
            path_info.append(root.graph.edges[(source, target)]['label'])
            if isinstance(self.elements[target], Regulation):
                path_info.append(self.elements[target].event.event_type)
            elif isinstance(self.elements[target], Entity):
                path_info.append(target)
        return path_info


def process_annotations(ann_file):
    """Iterate over the rows of an annotations file and build up objects."""
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


# Below we define dataclasses to represent elements of Standoff annotations

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
            v = [v] if not isinstance(v, list) else v
            for vv in v:
                if isinstance(vv, Regulation):
                    entrez_ids |= vv.find_entrez_ids()
                elif isinstance(vv, Entity):
                    entrez_id = vv.references.get('EG')
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
            paths.append(path_nodes)
        return paths


def add_subgraph(g, obj):
    """Recursively build up a graph of standoff objects."""
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
