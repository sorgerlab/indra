from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from indra.statements import *
from indra.databases import hgnc_client
import logging
import pybel
import pybel.constants as pc
from copy import deepcopy, copy
from pybel.parser.language import pmod_namespace
from pybel.parser.canonicalize import node_to_tuple
from indra.assemblers.pysb_assembler import mod_acttype_map

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('pybel_assembler')


_indra_pybel_act_map = {
    'kinase': 'kin',
    'phosphatase': 'phos',
    'catalytic': 'cat',
    'gtpbound': 'gtp',
    'transcription': 'tscript',
    'gef': 'gef',
    'gap': 'gap'
}


class PybelAssembler(object):
    def __init__(self, stmts=None, name=None, description=None, version=None,
                 **kwargs):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts

        # Create the model and assign metadata
        self.model = pybel.BELGraph(**kwargs)
        self.model.graph[pc.GRAPH_METADATA] = {
                pc.METADATA_NAME: name,
                pc.METADATA_DESCRIPTION: description,
                pc.METADATA_VERSION: version
            }

    def add_statements(self, stmts_to_add):
        self.stmts += stmts_to_add

    def make_model(self):
        for stmt in self.statements:
            # Skip statements with no subject
            if stmt.agent_list()[0] is None and \
                    not isinstance(stmt, Conversion):
                continue
            # Assemble statements
            if isinstance(stmt, Modification):
                self._assemble_modification(stmt)
            elif isinstance(stmt, RegulateActivity):
                self._assemble_regulate_activity(stmt)
            elif isinstance(stmt, RegulateAmount):
                self._assemble_regulate_amount(stmt)
            elif isinstance(stmt, Gef):
                self._assemble_gef(stmt)
            elif isinstance(stmt, Gap):
                self._assemble_gap(stmt)
            elif isinstance(stmt, ActiveForm):
                self._assemble_active_form(stmt)
            elif isinstance(stmt, Complex):
                self._assemble_complex(stmt)
            elif isinstance(stmt, Conversion):
                self._assemble_conversion(stmt)
            elif isinstance(stmt, Autophosphorylation):
                self._assemble_autophosphorylation(stmt)
            elif isinstance(stmt, Transphosphorylation):
                self._assemble_transphosphorylation(stmt)
            else:
                logger.info('Unhandled statement: %s' % stmt)
        return self.model

    def _assemble_regulate_activity(self, stmt):
        # Get node data and add to model
        subj_node, subj_attr, subj_edge = _get_agent_node(stmt.subj)
        act_obj = _get_activated_object(stmt)
        obj_node, obj_attr, obj_edge = _get_agent_node(act_obj)
        self.model.add_node(subj_node, attr_dict=subj_attr)
        self.model.add_node(obj_node, attr_dict=obj_attr)
        # Define the edge data
        pybel_relation = pc.DIRECTLY_INCREASES \
                         if isinstance(stmt, Activation) \
                         else pc.DIRECTLY_DECREASES
        edge_data_list = \
            _combine_edge_data(pybel_relation, subj_edge, obj_edge,
                               stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_modification(self, stmt):
        (_, enz_attr, enz_edge) = _get_agent_node(stmt.enz)
        sub_agent = _get_modified_substrate(stmt)
        (_, sub_attr, sub_edge) = _get_agent_node(sub_agent)
        enz_node = self.model.add_node_from_data(enz_attr)
        sub_node = self.model.add_node_from_data(sub_attr)
        pybel_relation = pc.DIRECTLY_INCREASES \
                         if isinstance(stmt, AddModification) \
                         else pc.DIRECTLY_DECREASES
        edge_data_list = _combine_edge_data(pybel_relation, enz_edge, sub_edge,
                                            stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(enz_node, sub_node, attr_dict=edge_data)

    def _assemble_regulate_amount(self, stmt):
        # p(HGNC:TP53) => p(HGNC:MDM2)
        (subj_node, subj_attr, subj_edge) = _get_agent_node(stmt.subj)
        (obj_node, obj_attr, obj_edge) = _get_agent_node(stmt.obj)
        self.model.add_node(subj_node, attr_dict=subj_attr)
        self.model.add_node(obj_node, attr_dict=obj_attr)
        pybel_relation = pc.DIRECTLY_INCREASES \
                         if isinstance(stmt, IncreaseAmount) \
                         else pc.DIRECTLY_DECREASES
        edge_data_list = _combine_edge_data(pybel_relation, subj_edge,
                                            obj_edge, stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_gef(self, stmt):
        gef = deepcopy(stmt.gef)
        gef.activity = ActivityCondition('gef', True)
        ras = deepcopy(stmt.ras)
        ras.activity = ActivityCondition('gtpbound', True)
        (subj_node, subj_attr, subj_edge) = _get_agent_node(gef)
        (obj_node, obj_attr, obj_edge) = _get_agent_node(ras)
        subj_node = self.model.add_node_from_data(subj_attr)
        obj_node = self.model.add_node_from_data(obj_attr)
        pybel_relation = pc.DIRECTLY_INCREASES
        edge_data_list = _combine_edge_data(pybel_relation, subj_edge,
                                            obj_edge, stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_gap(self, stmt):
        gap = deepcopy(stmt.gap)
        gap.activity = ActivityCondition('gap', True)
        ras = deepcopy(stmt.ras)
        ras.activity = ActivityCondition('gtpbound', True)
        (subj_node, subj_attr, subj_edge) = _get_agent_node(gap)
        (obj_node, obj_attr, obj_edge) = _get_agent_node(ras)
        self.model.add_node(subj_node, attr_dict=subj_attr)
        self.model.add_node(obj_node, attr_dict=obj_attr)
        pybel_relation = pc.DIRECTLY_DECREASES
        edge_data_list = _combine_edge_data(pybel_relation, subj_edge,
                                            obj_edge, stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_active_form(self, stmt):
        act_agent = Agent(stmt.agent.name, db_refs=stmt.agent.db_refs)
        act_agent.activity = ActivityCondition(stmt.activity, True)
        pybel_relation = pc.DIRECTLY_INCREASES if stmt.is_active \
                                               else pc.DIRECTLY_DECREASES
        (subj_node, subj_attr, subj_edge) = _get_agent_node(stmt.agent)
        (obj_node, obj_attr, obj_edge) = _get_agent_node(act_agent)
        self.model.add_node(subj_node, attr_dict=subj_attr)
        self.model.add_node(obj_node, attr_dict=obj_attr)
        edge_data_list = _combine_edge_data(pybel_relation, subj_edge,
                                            obj_edge, stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_complex(self, stmt):
        _, complex_node_data, _ = _get_complex_node(stmt.members)
        self.model.add_node_from_data(complex_node_data)

    def _assemble_conversion(self, stmt):
        pybel_lists = ([], [])
        for pybel_list, agent_list in \
                            zip(pybel_lists, (stmt.obj_from, stmt.obj_to)):
            for ag in agent_list:
                func, namespace, name = _get_agent_grounding(ag)
                pybel_list.append({
                    pc.FUNCTION: func,
                    pc.NAMESPACE: namespace,
                    pc.NAME: name})
        rxn_node_data = {
            pc.FUNCTION: pc.REACTION,
            pc.REACTANTS: pybel_lists[0],
            pc.PRODUCTS: pybel_lists[1]
        }
        obj_node = self.model.add_node_from_data(rxn_node_data)
        # Add node for controller, if there is one
        if stmt.subj is not None:
            # Add the node
            subj_node, subj_attr, subj_edge = _get_agent_node(stmt.subj)
            self.model.add_node_from_data(subj_attr)
            # Add the edge
            pybel_relation = pc.DIRECTLY_INCREASES
            edge_data_list = \
                _combine_edge_data(pybel_relation, subj_edge, None,
                                   stmt.evidence)
            for edge_data in edge_data_list:
                self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_autophosphorylation(self, stmt):
        (enz_node, enz_attr, enz_edge) = _get_agent_node(stmt.enz)
        sub_agent = deepcopy(stmt.enz)
        mc = stmt._get_mod_condition()
        sub_agent.mods.append(mc)
        (sub_node, sub_attr, sub_edge) = _get_agent_node(sub_agent)
        self.model.add_node(enz_node, attr_dict=enz_attr)
        self.model.add_node(sub_node, attr_dict=sub_attr)
        pybel_relation = pc.DIRECTLY_INCREASES
        edge_data_list = _combine_edge_data(pybel_relation, enz_edge, sub_edge,
                                            stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(enz_node, sub_node, attr_dict=edge_data)

    def _assemble_transphosphorylation(self, stmt):
        # Check our assumptions about the bound condition of the enzyme
        assert len(stmt.enz.bound_conditions) == 1
        assert stmt.enz.bound_conditions[0].is_bound
        # Get the "enzyme" node as the complex
        _, complex_data, _ = _get_agent_node(stmt.enz)
        # Create a modified protein node for the bound target
        sub_agent = deepcopy(stmt.enz.bound_conditions[0].agent)
        mc = stmt._get_mod_condition()
        sub_agent.mods.append(mc)
        (_, sub_attr, sub_edge) = _get_agent_node(sub_agent)
        complex_node = self.model.add_node_from_data(complex_data)
        sub_node = self.model.add_node_from_data(sub_attr)
        pybel_relation = pc.DIRECTLY_INCREASES
        edge_data_list = _combine_edge_data(pybel_relation, None, sub_edge,
                                            stmt.evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(complex_node, sub_node, attr_dict=edge_data)


def _combine_edge_data(relation, subj_edge, obj_edge, evidence):
    edge_data = {pc.RELATION: relation}
    if subj_edge:
        edge_data[pc.SUBJECT] = subj_edge
    if obj_edge:
        edge_data[pc.OBJECT] = obj_edge
    if not evidence:
        return [edge_data]
    edge_data_list = []
    for ev in evidence:
        pybel_ev = _get_evidence(ev)
        edge_data_one = copy(edge_data)
        edge_data_one.update(pybel_ev)
        edge_data_list.append(edge_data_one)
    return edge_data_list


def _get_modified_substrate(mod_stmt):
    mod_agent = deepcopy(mod_stmt.sub)
    mc = mod_stmt._get_mod_condition()
    mod_agent.mods.append(mc)
    return mod_agent


def _get_activated_object(reg_stmt):
    act_agent = deepcopy(reg_stmt.obj)
    ac = reg_stmt._get_activity_condition()
    # We set is_active to True here since the polarity is encoded
    # in the edge (decreases/increases)
    ac.is_active = True
    act_agent.activity = ac
    return act_agent


def _get_agent_node(agent):
    if agent.bound_conditions:
        members = [agent] + [bc.agent for bc in agent.bound_conditions
                             if bc.is_bound]
        return _get_complex_node(members)
    else:
        return _get_agent_node_no_bcs(agent)


def _get_complex_node(members):
    members_list = []
    for member in members:
        func, namespace, name = _get_agent_grounding(member)
        members_list.append({
            pc.FUNCTION: func,
            pc.NAMESPACE: namespace,
            pc.NAME: name})
    complex_node_data = {
            pc.FUNCTION: pc.COMPLEX,
            pc.MEMBERS: members_list}
    # TODO: Refactor to allow activity on a complex drawn from the prime agent
    # from a set of bound conditions
    return (None, complex_node_data, None)


def _get_agent_node_no_bcs(agent):
    (abundance_type, db_ns, db_id) = _get_agent_grounding(agent)
    if abundance_type is None:
        logger.warning('Agent %s has no grounding.', agent)
        return None
    node_attr = {pc.FUNCTION: abundance_type,
                 pc.NAMESPACE: db_ns,
                 pc.NAME: db_id}
    variants = []
    for mod in agent.mods:
        pybel_mod = pmod_namespace.get(mod.mod_type)
        if not pybel_mod:
            logger.info('Skipping modification of type %s on agent %s',
                         mod.mod_type, agent)
            continue
        var = {pc.KIND: pc.PMOD,
               pc.IDENTIFIER: {
                   pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE,
                   pc.NAME: pybel_mod}}
        if mod.residue is not None:
            res = amino_acids[mod.residue]['short_name'].capitalize()
            var[pc.PMOD_CODE] = res
        if mod.position is not None:
            var[pc.PMOD_POSITION] = int(mod.position)
        variants.append(var)
    for mut in agent.mutations:
        var = {pc.KIND: pc.HGVS,
               pc.IDENTIFIER: mut.to_hgvs()}
        variants.append(var)
    if variants:
        node_attr[pc.VARIANTS] = variants
    node_tuple = node_to_tuple(node_attr)
    # Also get edge data for the agent
    edge_data = _get_agent_activity(agent)
    return (node_tuple, node_attr, edge_data)


def _get_agent_grounding(agent):
    hgnc_id = agent.db_refs.get('HGNC')
    uniprot_id = agent.db_refs.get('UP')
    be_id = agent.db_refs.get('BE')
    pfam_id = agent.db_refs.get('PF')
    fa_id = agent.db_refs.get('FA')
    chebi_id = agent.db_refs.get('CHEBI')
    pubchem_id = agent.db_refs.get('PUBCHEM')
    go_id = agent.db_refs.get('GO')
    mesh_id = agent.db_refs.get('MESH')
    if hgnc_id:
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        return (pc.PROTEIN, 'HGNC', hgnc_name)
    elif uniprot_id:
        return (pc.PROTEIN, 'UP', uniprot_id)
    elif be_id:
        return (pc.PROTEIN, 'BE', be_id)
    elif pfam_id:
        return (pc.PROTEIN, 'PFAM', be_id)
    elif fa_id:
        return (pc.PROTEIN, 'NXPFA', be_id)
    elif chebi_id:
        return (pc.ABUNDANCE, 'CHEBI', chebi_id)
    elif pubchem_id:
        return (pc.ABUNDANCE, 'PUBCHEM', pubchem_id)
    elif go_id:
        return (pc.BIOPROCESS, 'GO', go_id)
    elif mesh_id:
        return (pc.BIOPROCESS, 'MESH', mesh_id)
    else:
        return (None, None, None)


def _get_agent_activity(agent):
    ac = agent.activity
    if not ac:
        return None
    if not ac.is_active:
        logger.warning('Cannot represent negative activity in PyBEL: %s' %
                       agent)
    edge_data = {pc.MODIFIER: pc.ACTIVITY}
    if not ac.activity_type == 'activity':
        pybel_activity = _indra_pybel_act_map[ac.activity_type]
        edge_data[pc.EFFECT] = {pc.NAME: pybel_activity,
                                pc.NAMESPACE: pc.BEL_DEFAULT_NAMESPACE}
    return edge_data


def _get_evidence(evidence):
    pybel_ev = {pc.EVIDENCE: evidence.text}
    # If there is a PMID, use it as the citation
    if evidence.pmid:
        citation = {pc.CITATION_TYPE: pc.CITATION_TYPE_PUBMED,
                    pc.CITATION_REFERENCE: evidence.pmid}
    # If no PMID, include the interface and source_api for now--
    # in general this should probably be in the annotations for all evidence
    else:
        cit_source = evidence.source_api if evidence.source_api else 'Unknown'
        cit_id = evidence.source_id if evidence.source_id else 'Unknown'
        cit_ref_str = '%s:%s' % (cit_source, cit_id)
        citation = {pc.CITATION_TYPE: pc.CITATION_TYPE_OTHER,
                    pc.CITATION_REFERENCE: cit_ref_str}
    pybel_ev[pc.CITATION] = citation
    pybel_ev[pc.ANNOTATIONS] = {}
    return pybel_ev


"""
Representation of PTM reactions in PyBEL/BEL
--------------------------------------------
- Kin -> pmod(Sub) vs.
- Kin -> rxn(Sub + ATP, pmod(Sub) + ADP)


structure of bel graph.
- node/tuple, maps to
  - dict with node/tuples as keys, each one mapped to
    - a dict, with integers (negative integers???) as keys, each one mapped to
      - a dict with edge information

inner dict has fields such as
{'annotations': {},
'citation': {
    'authors': '',
    'comments': '',
    'date': '',
    'name': 'Mol Cell Biol 2001 Apr 21(8) 2659-70',
    'reference': '11283246',
    'type': 'PubMed'},
'evidence': 'Modified assertion',
'subject': {
    'effect': {'name': 'kin', 'namespace': 'bel'},
    'modifier': 'Activity'}}}
'relation': 'increases',
'object': {
    'effect': {'name': 'kin', 'namespace': 'bel'},
    'modifier': 'Activity'},

Note that Activity modifiers are listed in the subject/object dicts.

Edges are added between top level nodes and modified nodes using a hasVariant
relation:

('Protein',
  'HGNC',
  'BRAF',
  ('hgvs', 'p.Val600Glu')): {-4: {'relation': 'hasVariant'}},


"""

# Need function for building up terms from agent conditions, including
# (mainly for now) modification conditions

# Example nodes in graph:
"""
('Protein', 'RGD', 'Pdpk1', ('pmod', ('bel', 'Ph'))),
 ('Complex', ('Protein', 'HGNC', 'JUP'), ('Protein', 'PFH', 'AXIN Family')),
 ('RNA', 'HGNC', 'LTV1'),
 ('RNA', 'MGI', 'Nfix'),
 ('Protein', 'HGNC', 'PAX7'),
 ('RNA', 'HGNC', 'NEUROG1'),
 ('RNA', 'MGI', 'Mup1'),
 ('Protein', 'HGNC', 'PPARA'),
 ('RNA', 'EGID', '66935'),
 ('RNA', 'HGNC', 'AKIRIN1'),
"""

# Then need patterns for assembling modifications, activations, rasgef/gap

# Can start with this, then add evidence/context, etc.


