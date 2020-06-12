import json
import logging
import itertools
import collections
import numpy as np
from copy import deepcopy
from indra.statements import *
from indra.databases import context_client, get_identifiers_url
from indra.tools.expand_families import Expander
from indra.ontology.bio import bio_ontology
from indra.ontology.standardize import \
    standardize_db_refs

expander = Expander(ontology=bio_ontology)


logger = logging.getLogger(__name__)


class CyJSAssembler(object):
    """This class assembles a CytoscapeJS graph from a set of INDRA Statements.

    CytoscapeJS is a web-based network library for analysis and
    visualisation: http://js.cytoscape.org/

    Parameters
    ----------
    statements : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    """
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
        self._gene_names = []
        self._context = {}

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
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

        Returns
        -------
        cyjs_str : str
            The json serialized Cytoscape JS model.
        """
        for stmt in self.statements:
            if isinstance(stmt, RegulateActivity):
                self._add_regulate_activity(stmt)
            elif isinstance(stmt, RegulateAmount):
                self._add_regulate_amount(stmt)
            elif isinstance(stmt, Modification):
                self._add_modification(stmt)
            elif isinstance(stmt, SelfModification):
                self._add_selfmodification(stmt)
            elif isinstance(stmt, Gef):
                self._add_gef(stmt)
            elif isinstance(stmt, Gap):
                self._add_gap(stmt)
            elif isinstance(stmt, Complex):
                self._add_complex(stmt)
            else:
                logger.warning('Unhandled statement type: %s' %
                               stmt.__class__.__name__)
        if kwargs.get('grouping'):
            self._group_nodes()
            self._group_edges()
        return self.print_cyjs_graph()

    def get_gene_names(self):
        """Gather gene names of all nodes and node members"""
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
        self._gene_names = gene_names

    def set_CCLE_context(self, cell_types):
        """Set context of all nodes and node members from CCLE."""
        self.get_gene_names()

        # Get expression and mutations from context client
        exp_values = \
            context_client.get_protein_expression(self._gene_names, cell_types)
        mut_values = \
            context_client.get_mutations(self._gene_names, cell_types)

        # Make a dict of presence/absence of mutations
        muts = {cell_line: {} for cell_line in cell_types}
        for cell_line, entries in mut_values.items():
            if entries is not None:
                for gene, mutations in entries.items():
                    if mutations:
                        muts[cell_line][gene] = 1
                    else:
                        muts[cell_line][gene] = 0

        # Create bins for the exp values
        # because colorbrewer only does 3-9 bins and I don't feel like
        # reinventing color scheme theory, this will only bin 3-9 bins
        def bin_exp(expression_dict):
            d = expression_dict
            exp_values = []
            for line in d:
                for gene in d[line]:
                    val = d[line][gene]
                    if val is not None:
                        exp_values.append(val)
            thr_dict = {}
            for n_bins in range(3, 10):
                bin_thr = np.histogram(np.log10(exp_values), n_bins)[1][1:]
                thr_dict[n_bins] = bin_thr
            # this dict isn't yet binned, that happens in the loop
            binned_dict = {x: deepcopy(expression_dict) for x in range(3, 10)}
            for n_bins in binned_dict:
                for line in binned_dict[n_bins]:
                    for gene in binned_dict[n_bins][line]:
                        # last bin is reserved for None
                        if binned_dict[n_bins][line][gene] is None:
                            binned_dict[n_bins][line][gene] = n_bins
                        else:
                            val = np.log10(binned_dict[n_bins][line][gene])
                            for thr_idx, thr in enumerate(thr_dict[n_bins]):
                                if val <= thr:
                                    binned_dict[n_bins][line][gene] = thr_idx
                                    break
            return binned_dict
        binned_exp = bin_exp(exp_values)

        context = {'bin_expression': binned_exp,
                   'mutation': muts}
        self._context['CCLE'] = context

    def print_cyjs_graph(self):
        """Return the assembled Cytoscape JS network as a json string.

        Returns
        -------
        cyjs_str : str
            A json string representation of the Cytoscape JS network.
        """
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        cyjs_str = json.dumps(cyjs_dict, indent=1, sort_keys=True)
        return cyjs_str

    def print_cyjs_context(self):
        """Return a list of node names and their respective context.

        Returns
        -------
        cyjs_str_context : str
            A json string of the context dictionary. e.g. -
            {'CCLE' : {'bin_expression' : {'cell_line1' : {'gene1':'val1'} },
            'bin_expression' : {'cell_line' : {'gene1':'val1'} }
            }}
        """
        context = self._context
        context_str = json.dumps(context, indent=1, sort_keys=True)
        return context_str

    def save_json(self, fname_prefix='model'):
        """Save the assembled Cytoscape JS network in a json file.

        This method saves two files based on the file name prefix given.
        It saves one json file with the graph itself, and another json
        file with the context.

        Parameters
        ----------
        fname_prefix : Optional[str]
            The prefix of the files to save the Cytoscape JS network and
            context to.
            Default: model
        """
        cyjs_str = self.print_cyjs_graph()
        # outputs the graph
        with open(fname_prefix + '.json', 'wb') as fh:
            fh.write(cyjs_str.encode('utf-8'))
        # outputs the context of graph nodes
        context_str = self.print_cyjs_context()
        with open(fname_prefix + '_context.json', 'wb') as fh:
            fh.write(context_str.encode('utf-8'))

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
        with open(fname, 'wb') as fh:
            fh.write(s.encode('utf-8'))

    def _add_binary_regulation(self, stmt):
        subj, obj = stmt.agent_list()
        if subj is None:
            return
        edge_type, edge_polarity = _get_stmt_type(stmt)
        source_id = self._add_node(subj, uuid=stmt.uuid)
        target_id = self._add_node(obj, uuid=stmt.uuid)
        self._add_edge(edge_type, source_id, target_id, edge_polarity,
                       stmt.uuid)

    _add_regulate_activity = _add_binary_regulation
    _add_regulate_amount = _add_binary_regulation
    _add_modification = _add_binary_regulation
    _add_gef = _add_binary_regulation
    _add_gap = _add_binary_regulation

    def _add_selfmodification(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        source_id = self._add_node(stmt.enz, uuid=stmt.uuid)
        self._add_edge(edge_type, source_id, source_id, edge_polarity,
                       stmt.uuid)

    def _add_complex(self, stmt):
        edge_type, edge_polarity = _get_stmt_type(stmt)
        for m1, m2 in itertools.combinations(stmt.members, 2):
            m1_id = self._add_node(m1, uuid=stmt.uuid)
            m2_id = self._add_node(m2, uuid=stmt.uuid)
            self._add_edge(edge_type, m1_id, m2_id, edge_polarity,
                           stmt.uuid)

    def _get_edge_dict(self):
        """Return a dict of edges.

        Keyed tuples of (i, source, target, polarity)
        with lists of edge ids [id1, id2, ...]
        """
        edge_dict = collections.defaultdict(lambda: [])
        if len(self._edges) > 0:
            for e in self._edges:
                data = e['data']
                key = tuple([data['i'], data['source'],
                            data['target'], data['polarity']])
                edge_dict[key] = data['id']
        return edge_dict

    def _add_edge(self, edge_type, source, target, edge_polarity, uuid):
        edge_dict = self._get_edge_dict()
        uuids = collections.defaultdict(lambda: [])
        edge = {'data': {'i': edge_type,
                         'source': source, 'target': target,
                         'polarity': edge_polarity}}
        data = edge['data']
        key = tuple([data['i'], data['source'],
                    data['target'], data['polarity']])
        if key in edge_dict:
            val = edge_dict[key]
            edge = [e for e in self._edges if e['data']['id'] == val][0]
        else:
            edge['data']['id'] = self._get_new_id()
            self._edges.append(edge)
        if type(uuid) is not list:
            uuid = [uuid]
        edge['data']['uuid_list'] = edge['data'].get('uuid_list', [])
        edge['data']['uuid_list'] += uuid
        return

    def _add_node(self, agent, uuid=None):
        node_key = agent.name
        node_id = self._existing_nodes.get(node_key)
        # if the node already exists we do not want to add it again
        # we must however add its uuid
        if node_id is not None:
            # fetch the appropriate node
            n = [x for x in self._nodes if x['data']['id'] == node_id][0]
            uuid_list = n['data']['uuid_list']
            if uuid not in uuid_list:
                uuid_list.append(uuid)
            return node_id
        db_refs = _get_db_refs(agent)
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node_name = agent.name
        node_name = node_name.replace('_', ' ')
        expanded_families = bio_ontology.get_children(*agent.get_grounding())
        expanded_families = [ch for ch in expanded_families if
                             ch[0] == 'HGNC']
        members = {}
        for member in expanded_families:
            member_db_refs = {member[0]: member[1]}
            member_db_refs = standardize_db_refs(member_db_refs)
            gene_name = bio_ontology.get_name(*member)
            members[gene_name] = {'db_refs': member_db_refs}
        node = {'data': {'id': node_id, 'name': node_name,
                         'db_refs': db_refs, 'parent': '',
                         'members': members, 'uuid_list': [uuid]}}
        self._nodes.append(node)
        return node_id

    def _get_new_id(self):
        ret = self._id_counter
        self._id_counter += 1
        return ret

    def _get_node_key(self, node_dict_item):
        """Return a tuple of sorted sources and targets given a node dict."""
        s = tuple(sorted(node_dict_item['sources']))
        t = tuple(sorted(node_dict_item['targets']))
        return (s, t)

    def _get_node_groups(self):
        """Return a list of node id lists that are topologically identical.

        First construct a node_dict which is keyed to the node id and
        has a value which is a dict with keys 'sources' and 'targets'.
        The 'sources' and 'targets' each contain a list of tuples
        (i, polarity, source) edge of the node. node_dict is then processed
        by _get_node_key() which returns a tuple of (s,t) where s,t are
        sorted tuples of the ids for the source and target nodes. (s,t) is
        then used as a key in node_key_dict where the values are the node
        ids. node_groups is restricted to groups greater than 1 node.
        """
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
        """Group all edges that are topologically identical.

        This means that (i, source, target, polarity) are the same, then sets
        edges on parent (i.e. - group) nodes to 'Virtual' and creates a new
        edge to represent all of them.
        """
        # edit edges on parent nodes and make new edges for them
        edges_to_add = [[], []]  # [group_edges, uuid_lists]
        for e in self._edges:
            new_edge = deepcopy(e)
            new_edge['data'].pop('id', None)
            uuid_list = new_edge['data'].pop('uuid_list', [])
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
            if source_node['data']['parent'] != '':
                new_edge['data']['source'] = source_node['data']['parent']
                e['data']['i'] = 'Virtual'
            # If the targete node is in a group, we change the target of this
            # edge to the group
            if target_node['data']['parent'] != '':
                new_edge['data']['target'] = target_node['data']['parent']
                e['data']['i'] = 'Virtual'
            if e['data']['i'] == 'Virtual':
                if new_edge not in edges_to_add[0]:
                    edges_to_add[0].append(new_edge)
                    edges_to_add[1].append(uuid_list)
                else:
                    idx = edges_to_add[0].index(new_edge)
                    edges_to_add[1][idx] += uuid_list
                    edges_to_add[1][idx] = list(set(edges_to_add[1][idx]))
        for ze in zip(*edges_to_add):
            edge = ze[0]
            edge['data']['id'] = self._get_new_id()
            edge['data']['uuid_list'] = ze[1]
            self._edges.append(edge)

    def _group_nodes(self):
        node_groups = self._get_node_groups()
        for group in node_groups:
            # Make new group node
            new_group_node = {'data': {'id': (self._get_new_id()),
                                       'name': ('Group' + str(group)),
                                       'parent': '', 'uuid_list': []}}
            member_nodes = [x for x in self._nodes if x['data']['id'] in group]
            for m_node in member_nodes:
                new_group_node['data']['uuid_list'] += \
                    m_node['data']['uuid_list']
                new_group_node['data']['uuid_list'] = \
                    list(set(new_group_node['data']['uuid_list']))
            # Point the node to its parent
            for node in self._nodes:
                if node['data']['id'] in group:
                    node['data']['parent'] = new_group_node['data']['id']
            self._nodes.append(new_group_node)


def _get_db_refs(agent):
    cyjs_db_refs = {}
    for db_name, db_ids in agent.db_refs.items():
        if isinstance(db_ids, int):
            db_id = str(db_ids)
        elif isinstance(db_ids, str):
            db_id = db_ids
        else:
            db_id = db_ids[0]
        if db_name == 'TEXT':
            url = db_id
        else:
            url = get_identifiers_url(db_name, db_id)
        if not url:
            continue
        db_name_map = {
            'UP': 'UniProt', 'PUBCHEM': 'PubChem',
            'IP': 'InterPro', 'NXPFA': 'NextProtFamily',
            'PF': 'Pfam', 'CHEBI': 'ChEBI'}
        name = db_name_map.get(db_name)
        if not name:
            name = db_name
        cyjs_db_refs[name] = url
    return cyjs_db_refs


def _get_stmt_type(stmt):
    if isinstance(stmt, AddModification):
        edge_type = stmt.__class__.__name__
        edge_polarity = 'positive'
    elif isinstance(stmt, RemoveModification):
        edge_type = stmt.__class__.__name__
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
    elif isinstance(stmt, DecreaseAmount):
        edge_type = 'DecreaseAmount'
        edge_polarity = 'negative'
    elif isinstance(stmt, IncreaseAmount):
        edge_type = 'IncreaseAmount'
        edge_polarity = 'positive'
    elif isinstance(stmt, Gef):
        edge_type = 'Gef'
        edge_polarity = 'positive'
    elif isinstance(stmt, Gap):
        edge_type = 'Gap'
        edge_polarity = 'negative'
    else:
        edge_type = stmt.__class__.__str__()
        edge_polarity = 'none'
    return edge_type, edge_polarity
