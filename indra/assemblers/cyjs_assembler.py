from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import json
import logging
import itertools
import collections
from indra.statements import *
from indra.databases import context_client
import indra.preassembler.hierarchy_manager as hm
from numpy import histogram
import indra.tools.expand_families as exp_fam
import indra.preassembler as pr

hierarchies = hm.hierarchies
expander = exp_fam.Expander(hierarchies)

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger('cyjs_assembler')

class CyJSAssembler(object):

    def __init__(self, stmts=None):
        if not stmts:
            self.statements = []
        else:
            self.statements = stmts
        self._edges = []
        self._nodes = []
        self._existing_nodes = {}
        self._id_counter = 0

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        stmts = pr.Preassembler.combine_duplicate_stmts(stmts)
        for stmt in stmts:
            self.statements.append(stmt)
    def make_model(self, *args, **kwargs):
        """Assemble a Cytoscape JS network from INDRA Statements.

        This method assembles a Cytoscape JS network from the set of INDRA
        Statements added to the assembler.

        Parameters
        ----------
        grouping : bool
            If True, the nodes with identical incoming and outgoing edges
            are grouped and the corresponding edges are merged.

        drop_virtual_edges : bool
            If True, the original edges which have been collected and made
            virtual are discarded. If these edges are discarded, they are
            not seen by the cytoscape.js layout algorithms.

        add_edge_weights : bool
            If True, give edges that connect group nodes a weight of their
            group size. All other edges get a weight of 1.

        Returns
        -------
        cyjs_str : str
            The json serialized Cytoscape JS model.
        """
        for stmt in self.statements:
            if isinstance(stmt, Activation):
                self._add_activation(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            elif isinstance(stmt, Modification):
                self._add_modification(stmt)
            else:
                logger.warning('Unhandled statement type: %s' %
                               stmt.__class__.__name__)
        if kwargs.get('grouping', None):
            self._group_nodes()
            self._group_edges()
        if kwargs.get('drop_virtual_edges', None):
            self._drop_virtual_edges()
        if kwargs.get('add_edge_weights', None):
            self._add_edge_weights()
        return self.print_cyjs()


    def set_context(self, *args, **kwargs):
        """Set protein expression data as node attribute

        This method uses :py:mod:`indra.databases.context_client` to get
        protein expression levels for a given cell type and set a node
        attribute for proteins accordingly.

        Parameters
        ----------
        cell_type : str
            Cell type name for which expression levels are queried.
            The cell type name follows the CCLE database conventions.

        Example: LOXIMVI_SKIN, BT20_BREAST
        """
        if kwargs.get('cell_type', None) :
            cell_type = kwargs.get('cell_type', None)
            node_names = [node['data']['name'] for node in self._nodes]
            exp = context_client.get_protein_expression(node_names, cell_type)
            mut = context_client.get_mutations(node_names, cell_type)
            if not exp:
                logger.warning('Could not get context for %s cell type.' %
                               cell_type)
                return
            if not mut:
                logger.warning('Could not get mutations for %s cell type.' %
                               cell_type)
                return
            counter_exp = 0
            counter_mut = 0
            for node in self._nodes:
                amount = exp.get(node['data']['name'])
                if amount is  None:
                    node['data']['expression'] = None
                if amount is not None:
                    node['data']['expression'] = int(amount[cell_type])
                    counter_exp += 1
                mutation = mut.get(node['data']['name'])
                if mutation is  None:
                    node['data']['mutation'] = None
                if mutation is not None:
                    node['data']['mutation'] = int(mutation[cell_type])
                    counter_mut += 1
            logger.info('Set expression context for %d nodes.' % counter_exp)
            logger.info('Set mutation context for %d nodes.' % counter_mut)
        if kwargs.get('bin_expression', None):
            exp_lvls = [n['data'].get('expression', None) for n in self._nodes]
            exp_lvls = [x for x in exp_lvls if x != None]
            bin_thr = histogram([x for x in exp_lvls if x != None], 5)[1]
            for n in self._nodes:
                if n['data']['expression'] ==  None:
                    n['data']['expression'] = 5
                else:
                    for thr_idx, thr in enumerate(bin_thr):
                        if n['data']['expression']<= thr:
                            n['data']['expression'] = thr_idx
                            break

    def print_cyjs(self):
        """Return the assembled Cytoscape JS network as a json string.

        Returns
        -------
        cyjs_str : str
            A json string representation of the Cytoscape JS network.
        """
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        cyjs_str = json.dumps(cyjs_dict, indent=1)
        return cyjs_str

    def save_model(self, fname='model.js'):
        """Save the assembled Cytoscape JS network in a file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the Cytoscape JS network to.
            Default: model.js
        """
        cyjs_str = self.print_cyjs()
        s = 'var modelElements = %s;' % cyjs_str
        with open(fname, 'wt') as fh:
            fh.write(s)

    def _add_activation(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        edge_id = self._get_new_id()
        source_id = self._add_node(stmt.subj)
        target_id = self._add_node(stmt.obj)
        edge = {'data': {'i': edge_type, 'id': edge_id,
                         'source': source_id, 'target': target_id,
                         'polarity': edge_polarity}}
        self._edges.append(edge)

    def _add_modification(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        edge_id = self._get_new_id()
        source_id = self._add_node(stmt.enz)
        target_id = self._add_node(stmt.sub)
        edge = {'data': {'i': edge_type, 'id': edge_id,
                         'source': source_id, 'target': target_id,
                         'polarity': edge_polarity}}
        self._edges.append(edge)

    def _add_complex(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self._add_node(m1)
            m2_id = self._add_node(m2)

            edge_id = self._get_new_id()
            edge = {'data': {'i': edge_type, 'id': edge_id,
                             'source': m1_id, 'target': m2_id,
                             'polarity': edge_polarity}}
            self._edges.append(edge)

    def _add_node(self, agent):
        node_key = agent.name
        node_id = self._existing_nodes.get(node_key)
        if node_id is not None:
            return node_id
        db_refs = _get_db_refs(agent)
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node_name = agent.name
        expanded_families = expander.get_children(agent, ns_filter='HGNC')
        members = {'HGNC':sorted([x[1] for x in expanded_families])}
        node = {'data': {'id': node_id, 'name': node_name,
                         'db_refs': db_refs, 'parent':'',
                         'members': members}}
        self._nodes.append(node)
        return node_id

    def _get_new_id(self):
        ret = self._id_counter
        self._id_counter += 1
        return ret

    def _get_node_key(self, node_dict):
        s = tuple(sorted(node_dict['sources']))
        t = tuple(sorted(node_dict['targets']))
        return (s, t)

    def _get_node_groups(self):
        # First we construct a dictionary for each node's
        # source and target edges
        node_dict = {node['data']['id']: {'sources': [], 'targets': []}
                     for node in self._nodes}
        for edge in self._edges:
            # Add edge as a source for its target node
            edge_data = (edge['data']['i'], edge['data']['polarity'],
                         edge['data']['source'])
            node_dict[edge['data']['target']]['sources'].append(edge_data)
            # Add edge as target for its source node
            edge_data = (edge['data']['i'], edge['data']['polarity'],
                         edge['data']['target'])
            node_dict[edge['data']['source']]['targets'].append(edge_data)

        # Make a dictionary of nodes based on source/target as a key
        node_key_dict = collections.defaultdict(lambda: [])
        for node_id, node_d in node_dict.items():
            key = self._get_node_key(node_d)
            node_key_dict[key].append(node_id)
        # Constrain the groups to ones that have more than 1 member
        node_groups = [g for g in node_key_dict.values() if (len(g) > 1)]
        return node_groups

    def _group_edges(self):
        # Iterate over edges in a copied edge list
        edges_to_add = []
        for e in self._edges:
            # Check if edge source or target are contained in a parent
            # If source or target in parent edit edge
            # Nodes may only point within their container
            source = e['data']['source']
            target = e['data']['target']
            source_node = [x for x in self._nodes if
                           x['data']['id'] == source][0]
            target_node = [x for x in self._nodes if
                           x['data']['id'] == target][0]
            # If the source node is in a group, we change the source of this
            # edge to the group
            new_edge = None
            if source_node['data']['parent'] != '':
                new_edge = deepcopy(e)
                new_edge['data'].pop('id', None)
                new_edge['data']['source'] = source_node['data']['parent']
                e['data']['i'] = 'Virtual'
            # If the targete node is in a group, we change the target of this
            # edge to the group
            if target_node['data']['parent'] != '':
                if new_edge is None:
                    new_edge = deepcopy(e)
                    new_edge['data'].pop('id', None)
                new_edge['data']['target'] = target_node['data']['parent']
                e['data']['i'] = 'Virtual'
            if new_edge is not None:
                if new_edge not in edges_to_add:
                    edges_to_add.append(new_edge)

        # need to check if there are identical edges in edges to add
        # identical on everything but id
        for edge in edges_to_add:
            new_id = self._get_new_id()
            edge['data']['id'] = new_id
            self._edges.append(edge)

    def _group_nodes(self):
        node_groups = self._get_node_groups()
        for group in node_groups:
            # Make new group node
            new_group_node = {'data': {'id': (self._get_new_id()),
                                       'name': ('Group' + str(group)),
                                       'parent': ''}}
            # Point the node to its parent
            for node in self._nodes:
                if node['data']['id'] in group:
                    node['data']['parent'] = new_group_node['data']['id']
            self._nodes.append(new_group_node)

    def _drop_virtual_edges(self):
        self._edges = [x for x in self._edges if x['data']['i'] != 'Virtual']

    def _add_edge_weights(self):
        # make a list of group nodes
        group_node_ids = []
        for n in self._nodes:
            if n['data']['parent'] != '':
                group_node_ids.append(n['data']['parent'])
        group_node_ids = list(set(group_node_ids))
        # get sizes for each group
        group_node_sizes = {}
        for g in group_node_ids:
            group_members = [x for x in self._nodes
                             if x['data']['parent'] == g]
            group_size = len(group_members)
            group_node_sizes[g] = group_size
        # iterate over edges
        # if they point to/from group, weigh them acc to group size
        # nodes between two groups get assigned heaviest of two weights
        for e in self._edges:
            source = e['data']['source']
            target = e['data']['target']
            if (source in group_node_ids) and (target in group_node_ids):
                e['data']['weight'] = max(group_node_sizes[source],
                                          group_node_sizes[target])
            elif source in group_node_ids:
                e['data']['weight'] = group_node_sizes[source]
            elif target in group_node_ids:
                e['data']['weight'] = group_node_sizes[target]
        # once all group node edges have weights
        # give non-group node edges weights of 1
        for e in self._edges:
            if e['data'].get('weight', None) is None:
                e['data']['weight'] = 1

def _get_db_refs(agent):
    cyjs_db_refs = {}
    for db_name, db_ids in agent.db_refs.items():
        if isinstance(db_ids, int):
            db_id = str(db_ids)
        elif isinstance(db_ids, basestring):
            db_id = db_ids
        else:
            db_id = db_ids[0]
        if db_name == 'UP':
            name = 'UniProt'
            val = 'http://identifiers.org/uniprot/%s' % db_id
        elif db_name == 'HGNC':
            name = 'HGNC'
            val = 'http://identifiers.org/hgnc/HGNC:%s' % db_id
        elif db_name == 'CHEBI':
            name = 'ChEBI'
            val = 'http://identifiers.org/chebi/%s' % db_id
        elif db_name == 'PUBCHEM':
            name = 'PubChem'
            val = 'http://identifiers.org/pubchem.compound/%s' % db_id
        elif db_name == 'HMDB':
            name = 'HMDB'
            val = 'http://identifiers.org/hmdb/%s' % db_id
        elif db_name == 'GO':
            name = 'GO'
            val = 'http://identifiers.org/go/%s' % db_id
        elif db_name == 'MESH':
            name = 'MESH'
            val = 'http://identifiers.org/mesh/%s' % db_id
        elif db_name == 'IP':
            name = 'InterPro'
            val = 'http://identifiers.org/interpro/%s' % db_id
        elif db_name == 'TEXT':
            continue
        else:
            val = db_id
            name = db_name
        cyjs_db_refs[name] = val
    return cyjs_db_refs

def _get_stmt_type(stmt):
    if isinstance(stmt, Modification):
        edge_type = 'Modification'
        edge_polarity = 'positive'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 'positive'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 'none'
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
        if stmt.is_activation:
            edge_polarity = 'positive'
        else:
            edge_polarity = 'negative'
    elif isinstance(stmt, RasGef):
        edge_type = 'RasGef'
        edge_polarity = 'positive'
    elif isinstance(stmt, RasGap):
        edge_type = 'RasGap'
        edge_polarity = 'negative'
    else:
        edge_type = stmt.__class__.__str__()
        edge_polarity = 'none'
    return edge_type, edge_polarity
