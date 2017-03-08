from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from copy import deepcopy
import json
import logging
import itertools
import collections
import numpy as np
from matplotlib.colors import LinearSegmentedColormap as colormap
from matplotlib.colors import rgb2hex, hex2color
from indra.statements import *
from indra.databases import hgnc_client
from indra.databases import context_client
from indra.preassembler import Preassembler
from indra.tools.expand_families import Expander
from indra.preassembler.hierarchy_manager import hierarchies

expander = Expander(hierarchies)

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
        self._exp_colorscale = []
        self._mut_colorscale = []

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        stmts = Preassembler.combine_duplicate_stmts(stmts)
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
            if isinstance(stmt, RegulateActivity):
                self._add_regulate_activity(stmt)
            elif isinstance(stmt, Inhibition):
                self._add_activation(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            elif isinstance(stmt, Modification):
                self._add_modification(stmt)
            else:
                logger.warning('Unhandled statement type: %s' %
                               stmt.__class__.__name__)
        if kwargs.get('grouping'):
            self._group_nodes()
            self._group_edges()
        if kwargs.get('drop_virtual_edges'):
            self._drop_virtual_edges()
        if kwargs.get('add_edge_weights'):
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

        bin_expression : bool
            If True, the gene expression will be put into 5 bins based on
            all gene expression values. An additional bin is used to indicate
            that the context_client returned None.

        user_bins : int
            If specified, split the expression levels into the given number
            of bins. If not specified, default will be 5.
        """
        cell_type = kwargs.get('cell_type')
        if not cell_type:
            logger.warning('No cell type given.')
            return

        # Collect all gene names in network
        gene_names = []
        for node in self._nodes:
            members = node['data'].get('members')
            if members:
                gene_names += list(members.keys())
            else:
                if node['data']['name'].startswith('Group'):
                    continue
                gene_names.append(node['data']['name'])

        # Get expression and mutation from context client
        exp = context_client.get_protein_expression(gene_names, cell_type)
        mut = context_client.get_mutations(gene_names, cell_type)
        if not exp:
            logger.warning('Could not get context for %s cell type.' %
                           cell_type)
            return
        else:
            exp = {k: v[cell_type] for k, v in exp.items()}
        if not mut:
            logger.warning('Could not get mutations for %s cell type.' %
                           cell_type)
            return
        else:
            mut = {k: v[cell_type] for k, v in mut.items()}

        # Get expression and mutation for specific gene
        def get_expr_mut(name, expr_data, mut_data):
            amount = expr_data.get(name)
            if amount is None:
                expression = None
            else:
                expression = np.log10(amount)
            mutation = mut_data.get(name)
            if mutation is not None:
                mutation = int(mutation)
            else:
                mutation = 0
            return expression, mutation

        # Set node properties for expression and mutation
        for node in self._nodes:
            members = node['data'].get('members')
            if members:
                for member in members.keys():
                    expression, mutation = get_expr_mut(member, exp, mut)
                    node['data']['members'][member]['expression'] = expression
                    node['data']['members'][member]['mutation'] = mutation
                node['data']['expression'] = None
                node['data']['mutation'] = 0
            else:
                if node['data']['name'].startswith('Group'):
                    node['data']['expression'] = None
                    node['data']['mutation'] = 0
                else:
                    expression, mutation = get_expr_mut(node['data']['name'],
                                                        exp, mut)
                    node['data']['expression'] = expression
                    node['data']['mutation'] = mutation

        # Binning for the purpose of assigning colors
        if kwargs.get('bin_expression'):
            # how many bins? If not specified, set to 5
            n_bins = 5
            user_bins = kwargs.get('n_bins')
            if type(user_bins) == int:
                n_bins = user_bins
                if n_bins > 9:
                    n_bins = 9
                    logger.info('Only 9 bins allowed. Setting n_bins = 9.')
                if n_bins < 3:
                    n_bins = 3
                    logger.info('Need at least 3 bin. Setting n_bins = 3.')
            # Create color scale for unmutated gene expression
            # feed in hex values from colorbrewer2 9-class PuBuGn
            wt_hexes = ['#f7fcf5', '#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476',
                        '#41ab5d', '#238b45', '#006d2c', '#00441b']
            exp_wt_colorscale = _build_color_scale(wt_hexes, n_bins)
            # tack on a gray for no expression data
            exp_wt_colorscale.append('#bdbdbd')
            self._exp_colorscale = exp_wt_colorscale
            # create color scale for mutated gene expression
            # feed in hex values from colorbrewer2 9-class YlOrRd
            mut_hexes = ['#fff5eb', '#fee6ce', '#fdd0a2', '#fdae6b', '#fd8d3c',
                         '#f16913', '#d94801', '#a63603', '#7f2704']
            exp_mut_colorscale = _build_color_scale(mut_hexes, n_bins)
            # tack on a gray for no expression data
            exp_mut_colorscale.append('#bdbdbd')
            self._mut_colorscale = exp_mut_colorscale
            # capture the expression levels of every gene in nodes
            exp_lvls = [n['data'].get('expression') for n in self._nodes]
            # capture the expression levels of every gene in family members
            m_exp_lvls = []
            for n in self._nodes:
                if n['data'].get('members'):
                    members = n['data']['members']
                    for m in members:
                        m_exp_lvls.append(members[m]['expression'])
            # combine node expressions and family expressions
            exp_lvls = exp_lvls + m_exp_lvls
            # get rid of None gene expressions
            exp_lvls = [x for x in exp_lvls if x is not None]
            # bin expression levels into n equally sized bins
            # bin n+1 reserved for None
            # this returns the bounds of each bin. so n_bins+1 bounds.
            # get rid of first value which is the leftmost bound
            bin_thr = np.histogram(exp_lvls, n_bins)[1][1:]
            # iterate over nodes
            for n in self._nodes:
                # if node has members set member bin_expression values
                if n['data'].get('members'):
                    members = n['data']['members']
                    for m in members:
                        # if expression is None, set to bin index n_bins
                        # if there is an expression value, bin and set to None
                        if members[m]['expression'] is None:
                            members[m]['bin_expression'] = n_bins
                        else:
                            for thr_idx, thr in enumerate(bin_thr):
                                if members[m]['expression'] <= thr:
                                    members[m]['bin_expression'] = thr_idx
                                    members[m]['expression'] = None
                                    break
                # set bin_expression for the node itself
                # if there is an expression value, bin and set to None
                if n['data']['expression'] is None:
                    n['data']['bin_expression'] = n_bins
                else:
                    for thr_idx, thr in enumerate(bin_thr):
                        if n['data']['expression'] <= thr:
                            n['data']['bin_expression'] = thr_idx
                            n['data']['expression'] = None
                            break

    def print_cyjs(self):
        """Return the assembled Cytoscape JS network as a json string.

            Returns
            -------
            cyjs_str : str
            A json string representation of the Cytoscape JS network.
        """
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        model_dict = {'exp_colorscale': self._exp_colorscale,
                      'mut_colorscale': self._mut_colorscale,
                      'model_elements': cyjs_dict}
        cyjs_str = json.dumps(model_dict, indent=1, sort_keys=True)
        return cyjs_str

    def save_json(self, fname='model.json'):
        """Save the assembled Cytoscape JS network in a json file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the Cytoscape JS network to.
            Default: model.json
        """
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        model_dict = {'exp_colorscale': self._exp_colorscale,
                      'mut_colorscale': self._mut_colorscale,
                      'model_elements': cyjs_dict}
        json_str = json.dumps(model_dict, indent=1, sort_keys=True)
        with open(fname, 'wt') as fh:
            fh.write(json_str)

    def save_model(self, fname='model.js'):
        """Save the assembled Cytoscape JS network in a js file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the Cytoscape JS network to.
            Default: model.js
        """
        exp_colorscale_str = json.dumps(self._exp_colorscale)
        mut_colorscale_str = json.dumps(self._mut_colorscale)
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        model_str = json.dumps(cyjs_dict, indent=1, sort_keys=True)
        model_dict = {'exp_colorscale_str': exp_colorscale_str,
                      'mut_colorscale_str': mut_colorscale_str,
                      'model_elements_str': model_str}
        s = ''
        s += 'var exp_colorscale = %s;\n' % model_dict['exp_colorscale_str']
        s += 'var mut_colorscale = %s;\n' % model_dict['mut_colorscale_str']
        s += 'var model_elements = %s;\n' % model_dict['model_elements_str']
        with open(fname, 'wt') as fh:
            fh.write(s)

    def _add_regulate_activity(self, stmt):
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
        node_name = node_name.replace('_', ' ')
        expanded_families = expander.get_children(agent, ns_filter='HGNC')
        members = {}
        for member in expanded_families:
            hgnc_symbol = member[1]
            hgnc_id = hgnc_client.get_hgnc_id(hgnc_symbol)
            if hgnc_id:
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                member_agent = Agent(hgnc_symbol,
                                     db_refs={'HGNC': hgnc_id,
                                              'UP': up_id})
                member_db_refs = _get_db_refs(member_agent)
            else:
                member_db_refs = {}
            members[member[1]] = {
                    'mutation': None,
                    'expression': None,
                    'db_refs': member_db_refs
                    }
        node = {'data': {'id': node_id, 'name': node_name,
                         'db_refs': db_refs, 'parent': '',
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
            if e['data'].get('weight') is None:
                e['data']['weight'] = 1

    def _add_attractors(self):
        parent_node_ids = [x['data']['parent'] for x in self._nodes
                           if x['data']['parent'] != '']
        parent_node_ids = list(set(parent_node_ids))
        attr_dict = {}
        for parent_node_id in parent_node_ids:
            child_node_ids = [x['data']['id'] for x in self._nodes
                              if x['data']['parent'] == parent_node_id]
            #  pick the middle node of the children
            #  this actually produces some spread group nodes, not ideal
            for i in list(range(len(child_node_ids))):
                if i >= len(child_node_ids)/2:
                    break
                else:
                    attr_node_id = child_node_ids[i]
            # sets attractor to last child node
            # attr_node_id = child_node_ids[-1]
            attr_dict[parent_node_id] = attr_node_id
        # for any existing edges to/from parent
        # give attractors same edges
        attr_edges = []
        for edge in self._edges:
            source = edge['data']['source']
            target = edge['data']['target']
            attr_edge = None
            # check source and target against attr_dict to point to attractors
            # edges sourcing or targeting parents will be ignored
            if source in attr_dict or target in attr_dict:
                attr_edge = deepcopy(edge)
                if source in attr_dict:
                    attr_edge['data']['source'] = attr_dict[source]
                if target in attr_dict:
                    attr_edge['data']['target'] = attr_dict[target]
                attr_edge['data']['id'] = self._get_new_id()
                attr_edge['data']['i'] = 'Attractor'
            if attr_edge is not None:
                if attr_edge not in attr_edges:
                    attr_edges.append(attr_edge)
        for attr_edge in attr_edges:
            self._edges.append(attr_edge)

    def _add_ext_attractors(self):
        parent_node_ids = [x['data']['parent'] for x in self._nodes
                           if x['data']['parent'] != '']
        parent_node_ids = list(set(parent_node_ids))
        # this is in the format {parent_node_id : [child_node_ids]}
        # parent child dict
        pc_dict = {}
        for parent_node_id in parent_node_ids:
            child_node_ids = [x['data']['id'] for x in self._nodes
                              if x['data']['parent'] == parent_node_id]
            pc_dict[parent_node_id] = {'children': child_node_ids,
                                       'sources': [],
                                       'src_attr_id': None,
                                       'targets': [],
                                       'targ_attr_id': None}
        # discover all sources and targets for group nodes
        for e in self._edges:
            source = e['data']['source']
            target = e['data']['target']
            if source in pc_dict or target in pc_dict:
                # any edge that has a parent node as its source is a target
                # for that parent node
                if source in pc_dict:
                    pc_dict[source]['targets'].append(target)
                # any edge that has a parent node as a target is a source
                # for that parent node
                if target in pc_dict:
                    pc_dict[target]['sources'].append(source)
        # create external attractor nodes for each parent node
        for p in pc_dict:
            # if there are sources that point at the parent node
            # init and append a source attractor
            children = pc_dict[p]['children']
            sources = pc_dict[p]['sources']
            if len(sources) > 0:
                src_attr_id = self._get_new_id()
                pc_dict[p]['srt_attr_id'] = src_attr_id
                src_attr = {'data': {'id': src_attr_id,
                                     'name': ('Attractor'),
                                     'parent': ''}}
                self._nodes.append(src_attr)
                # create edges from the sources to the source attractor
                for s in sources:
                    edge = {'data': {'i': 'Attractor',
                                     'id': self._get_new_id(),
                                     'source': s,
                                     'target': src_attr_id}}
                    self._edges.append(edge)
                # create edges from the src attractor pointing to children
                for c in children:
                    edge = {'data': {'i': 'Attractor',
                                     'id': self._get_new_id(),
                                     'source': src_attr_id,
                                     'target': c}}
                    self._edges.append(edge)
            # if there are nodes targeted by the parent node
            # init and append a target attractor
            targets = pc_dict[p]['targets']
            if len(targets) > 0:
                targ_attr_id = self._get_new_id()
                pc_dict[p]['targ_attr_id'] = src_attr_id
                targ_attr = {'data': {'id': targ_attr_id,
                                      'name': ('Attractor'),
                                      'parent': ''}}
                self._nodes.append(targ_attr)
                # create edges from the target attractor to targets
                for t in targets:
                    edge = {'data': {'i': 'Attractor',
                                     'id': self._get_new_id(),
                                     'source': targ_attr_id,
                                     'target': t}}
                    self._edges.append(edge)
                # create edges from the src attractor pointing to children
                for c in children:
                    edge = {'data': {'i': 'Attractor',
                                     'id': self._get_new_id(),
                                     'source': c,
                                     'target': targ_attr_id}}
                    self._edges.append(edge)


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
    if isinstance(stmt, AddModification):
        edge_type = 'Modification'
        edge_polarity = 'positive'
    elif isinstance(stmt, RemoveModification):
        edge_type = 'Modification'
        edge_polarity = 'negative'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
        edge_polarity = 'positive'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
        edge_polarity = 'none'
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
        edge_polarity = 'positive'
    elif isinstance(stmt, Inhibition):
        edge_type = 'Inhibition'
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


def _build_color_scale(hex_colors_list, n_bins):
    rgb_colors = [hex2color(x) for x in hex_colors_list]
    rgb_colors_array = np.array(rgb_colors)
    rgb_names = {'red': 0, 'green': 1, 'blue': 2}
    linear_mapping = np.linspace(0, 1, len(rgb_colors_array))
    cdict = {}
    for rgb_name in rgb_names:
        color_list = []
        rgb_idx = rgb_names[rgb_name]
        for lin, val in zip(linear_mapping, rgb_colors_array[:, rgb_idx]):
            color_list.append((lin, val, val))
        cdict[rgb_name] = color_list
    cmap = colormap('expression_colormap', cdict, 256, 1)
    color_scale = []
    for i in np.linspace(0, 1, n_bins):
        color_scale.append(rgb2hex(cmap(i)))

    return color_scale
