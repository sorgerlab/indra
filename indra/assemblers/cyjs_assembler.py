from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import json
from indra.statements import *

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
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self):
        for stmt in self.statements:
            if isinstance(stmt, Activation):
                self._add_activation(stmt)

    def print_cyjs(self):
        cyjs_dict = {'edges': self._edges, 'nodes': self._nodes}
        cyjs_str = json.dumps(cyjs_dict, indent=1)
        return cyjs_str

    def _add_activation(self, stmt):
        edge_type = _get_stmt_type(stmt)
        edge_id = self._get_new_id()
        source_id = self._add_node(stmt.subj)
        target_id = self._add_node(stmt.obj)
        edge = {'data': {'i': edge_type, 'id': edge_id,
                         'source': source_id, 'target': target_id}}
        self._edges.append(edge)

    def _add_node(self, agent):
        node_key = agent.name
        node_id = self._existing_nodes.get(node_key)
        if node_id is not None:
            return node_id
        node_id = self._get_new_id()
        self._existing_nodes[node_key] = node_id
        node_name = agent.name
        node = {'data': {'id': node_id, 'name': node_name}}
        self._nodes.append(node)
        return node_id

    def _get_new_id(self):
        ret = self._id_counter
        self._id_counter += 1
        return ret

def _get_stmt_type(stmt):
    if isinstance(stmt, Modification):
        edge_type = 'Modification'
    elif isinstance(stmt, SelfModification):
        edge_type = 'SelfModification'
    elif isinstance(stmt, Complex):
        edge_type = 'Complex'
    elif isinstance(stmt, Activation):
        edge_type = 'Activation'
    elif isinstance(stmt, RasGef):
        edge_type = 'RasGef'
    elif isinstance(stmt, RasGap):
        edge_type = 'RasGap'
    else:
        edge_type = stmt.__class__.__str__()
    return edge_type
