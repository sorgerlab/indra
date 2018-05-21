from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import uuid
import logging
import networkx as nx
from copy import deepcopy, copy
import pybel
import pybel.constants as pc
from pybel.language import pmod_namespace
from indra.statements import *
from indra.databases import hgnc_client

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


_pybel_indra_act_map = {v: k for k, v in _indra_pybel_act_map.items()}


class PybelAssembler(object):
    """Assembles a PyBEL graph from a set of INDRA Statements.

    PyBEL tools can subsequently be used to export the PyBEL graph into BEL
    script files, SIF files, and other related output formats.

    Parameters
    ----------
    stmts : list[:py:class:`indra.statement.Statement`]
        The list of Statements to assemble.
    name : str
        Name of the assembled PyBEL network.
    description : str
        Description of the assembled PyBEL network.
    version : str
        Version of the assembled PyBEL network.
    authors : str
        Author(s) of the network.
    contact : str
        Contact information (email) of the responsible author.
    license : str
        License information for the network.
    copyright : str
        Copyright information for the network.
    disclaimer : str
        Any disclaimers for the network.

    Examples
    --------
    >>> from indra.statements import *
    >>> map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    >>> mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    >>> stmt = Phosphorylation(map2k1, mapk1, 'T', '185')
    >>> pba = PybelAssembler([stmt])
    >>> belgraph = pba.make_model()
    >>> sorted(belgraph.nodes()) # doctest:+IGNORE_UNICODE
    [('Protein', 'HGNC', 'MAP2K1'), ('Protein', 'HGNC', 'MAPK1'), ('Protein', 'HGNC', 'MAPK1', ('pmod', ('bel', 'Ph'), 'Thr', 185))]
    >>> len(belgraph)
    3
    >>> belgraph.number_of_edges()
    2
    """
    def __init__(self, stmts=None, name=None, description=None, version=None,
                 authors=None, contact=None, license=None, copyright=None,
                 disclaimer=None):
        if stmts is None:
            self.statements = []
        else:
            self.statements = stmts
        if name is None:
            name = 'indra'
        if version is None:
            version = str(uuid.uuid4())
        # Create the model and assign metadata
        self.model = pybel.BELGraph(
            name=name,
            description=description,
            version=version,
            authors=authors,
            contact=contact,
            license=license,
            copyright=copyright,
            disclaimer=disclaimer,
        )
        ns_dict = {
            'HGNC': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                    'namespace/hgnc-human-genes/hgnc-human-genes-20170725.belns',
            'UP': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                  'namespace/swissprot/swissprot-20170725.belns',
            'IP': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                  'namespace/interpro/interpro-20170731.belns',
            #'FPLX':
            #'PFAM':
            #'NXPFA':
            'CHEBI': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                     'namespace/chebi-ids/chebi-ids-20170725.belns',
            'GO': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                  'namespace/go/go-20180109.belns',
            'MESH': 'https://arty.scai.fraunhofer.de/artifactory/bel/'
                    'namespace/mesh-processes/mesh-processes-20170725.belns'
        }
        self.model.namespace_url.update(ns_dict)
        self.model.namespace_pattern['PUBCHEM'] = '\d+'

    def add_statements(self, stmts_to_add):
        self.statements += stmts_to_add

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

    def to_database(self, connection=None):
        """Send the model to the PyBEL database

        This function wraps :py:func:`pybel.to_database`.

        Parameters
        ----------
        connection : Optional[str or pybel.manager.Manager]
            An RFC-1738 database connection string to the PyBEL SQL database
            or a PyBEL manager. If none, first checks the PyBEL configuration
            for ``PYBEL_CONNECTION`` then checks the environment variable
            ``PYBEL_REMOTE_HOST``. Finally, defaults to using SQLite
            database in PyBEL data directory (automatically configured
            by PyBEL)

        Returns
        -------
        network : Optional[pybel.manager.models.Network]
            The SQLAlchemy model representing the network that was uploaded.
            Returns None if upload fails.
        """
        network = pybel.to_database(self.model, connection=connection)
        return network

    def to_web(self, host=None, user=None, password=None):
        """Send the model to BEL Commons by wrapping :py:func:`pybel.to_web`

        The parameters ``host``, ``user``, and ``password`` all check the
        PyBEL configuration, which is located at
        ``~/.config/pybel/config.json`` by default

        Parameters
        ----------
        host : Optional[str]
            The host name to use. If none, first checks the PyBEL
            configuration entry ``PYBEL_REMOTE_HOST``, then the
            environment variable ``PYBEL_REMOTE_HOST``. Finally, defaults
            to https://bel-commons.scai.fraunhofer.de.
        user : Optional[str]
            The username (email) to use. If none, first checks the
            PyBEL configuration entry ``PYBEL_REMOTE_USER``,
            then the environment variable ``PYBEL_REMOTE_USER``.
        password : Optional[str]
            The password to use. If none, first checks the PyBEL configuration
            entry ``PYBEL_REMOTE_PASSWORD``, then the environment variable
            ``PYBEL_REMOTE_PASSWORD``.

        Returns
        -------
        response : requests.Response
            The response from the BEL Commons network upload endpoint.
        """
        response = pybel.to_web(self.model, host=host, user=user,
                                password=password)
        return response

    def save_model(self, path, output_format=None):
        """Save the :class:`pybel.BELGraph` using one of the outputs from
        :py:mod:`pybel`

        Parameters
        ----------
        path : str
            The path to output to
        output_format : Optional[str]
            Output format as ``cx``, ``pickle``, ``json`` or defaults to ``bel``
        """
        if output_format == 'pickle':
            pybel.to_pickle(self.model, path)
        else:
            with open(path, 'w') as fh:
                if output_format == 'json':
                    pybel.to_json_file(self.model, fh)
                elif output_format == 'cx':
                    pybel.to_cx_file(self.model, fh)
                else: # output_format == 'bel':
                    pybel.to_bel(self.model, fh)


    def to_signed_graph(self, symmetric_variant_links=False):
        edge_set = set()
        for u, v, edge_data in self.model.edges(data=True):
            rel = edge_data.get('relation')
            if rel in (pc.INCREASES, pc.DIRECTLY_INCREASES):
                edge_set.add((u, v, 0))
            elif rel in (pc.HAS_VARIANT, pc.HAS_COMPONENT):
                edge_set.add((u, v, 0))
                if symmetric_variant_links:
                    edge_set.add((v, u, 0))
            elif rel in (pc.DECREASES, pc.DIRECTLY_DECREASES):
                edge_set.add((u, v, 1))
            else:
                continue
        # Turn the tuples into dicts
        edge_data = [(u, v, dict([('sign', sign)])) for u, v, sign in edge_set]
        graph = nx.MultiDiGraph()
        graph.add_edges_from(edge_data)
        return graph


    def _add_nodes_edges(self, subj_agent, obj_agent, relation, evidence):
        """Given subj/obj agents, relation, and evidence, add nodes/edges."""
        subj_data, subj_edge = _get_agent_node(subj_agent)
        obj_data, obj_edge = _get_agent_node(obj_agent)
        # If we failed to create nodes for subject or object, skip it
        if subj_data is None or obj_data is None:
            return
        subj_node = self.model.add_node_from_data(subj_data)
        obj_node = self.model.add_node_from_data(obj_data)
        edge_data_list = \
            _combine_edge_data(relation, subj_edge, obj_edge, evidence)
        for edge_data in edge_data_list:
            self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_regulate_activity(self, stmt):
        """Example: p(HGNC:MAP2K1) => act(p(HGNC:MAPK1))"""
        act_obj = deepcopy(stmt.obj)
        act_obj.activity = stmt._get_activity_condition()
        # We set is_active to True here since the polarity is encoded
        # in the edge (decreases/increases)
        act_obj.activity.is_active = True
        relation = pc.DIRECTLY_INCREASES if isinstance(stmt, Activation) \
                                         else pc.DIRECTLY_DECREASES
        self._add_nodes_edges(stmt.subj, act_obj, relation, stmt.evidence)

    def _assemble_modification(self, stmt):
        """Example: p(HGNC:MAP2K1) => p(HGNC:MAPK1, pmod(Ph, Thr, 185))"""
        sub_agent = deepcopy(stmt.sub)
        sub_agent.mods.append(stmt._get_mod_condition())
        relation = pc.DIRECTLY_INCREASES if isinstance(stmt, AddModification) \
                                         else pc.DIRECTLY_DECREASES
        self._add_nodes_edges(stmt.enz, sub_agent, relation, stmt.evidence)

    def _assemble_regulate_amount(self, stmt):
        """Example: p(HGNC:ELK1) => p(HGNC:FOS)"""
        relation = pc.DIRECTLY_INCREASES if isinstance(stmt, IncreaseAmount) \
                                         else pc.DIRECTLY_DECREASES
        self._add_nodes_edges(stmt.subj, stmt.obj, relation, stmt.evidence)

    def _assemble_gef(self, stmt):
        """Example: act(p(HGNC:SOS1), ma(gef)) => act(p(HGNC:KRAS), ma(gtp))"""
        gef = deepcopy(stmt.gef)
        gef.activity = ActivityCondition('gef', True)
        ras = deepcopy(stmt.ras)
        ras.activity = ActivityCondition('gtpbound', True)
        self._add_nodes_edges(gef, ras, pc.DIRECTLY_INCREASES, stmt.evidence)

    def _assemble_gap(self, stmt):
        """Example: act(p(HGNC:RASA1), ma(gap)) =| act(p(HGNC:KRAS), ma(gtp))"""
        gap = deepcopy(stmt.gap)
        gap.activity = ActivityCondition('gap', True)
        ras = deepcopy(stmt.ras)
        ras.activity = ActivityCondition('gtpbound', True)
        self._add_nodes_edges(gap, ras, pc.DIRECTLY_DECREASES, stmt.evidence)

    def _assemble_active_form(self, stmt):
        """Example: p(HGNC:ELK1, pmod(Ph)) => act(p(HGNC:ELK1), ma(tscript))"""
        act_agent = Agent(stmt.agent.name, db_refs=stmt.agent.db_refs)
        act_agent.activity = ActivityCondition(stmt.activity, True)
        relation = pc.DIRECTLY_INCREASES if stmt.is_active \
                                         else pc.DIRECTLY_DECREASES
        self._add_nodes_edges(stmt.agent, act_agent, relation, stmt.evidence)

    def _assemble_complex(self, stmt):
        """Example: complex(p(HGNC:MAPK14), p(HGNC:TAB1))"""
        complex_data, _ = _get_complex_node(stmt.members)
        if complex_data is None:
            logger.info('skip adding complex with no members: %s', stmt.members)
            return
        self.model.add_node_from_data(complex_data)

    def _assemble_conversion(self, stmt):
        """Example: p(HGNC:HK1) => rxn(reactants(a(CHEBI:"CHEBI:17634")),
                                       products(a(CHEBI:"CHEBI:4170")))"""
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
        obj_edge = None # TODO: Any edge information possible here?
        # Add node for controller, if there is one
        if stmt.subj is not None:
            subj_attr, subj_edge = _get_agent_node(stmt.subj)
            subj_node = self.model.add_node_from_data(subj_attr)
            edge_data_list = _combine_edge_data(pc.DIRECTLY_INCREASES,
                                           subj_edge, obj_edge, stmt.evidence)
            for edge_data in edge_data_list:
                self.model.add_edge(subj_node, obj_node, attr_dict=edge_data)

    def _assemble_autophosphorylation(self, stmt):
        """Example: complex(p(HGNC:MAPK14), p(HGNC:TAB1)) =>
                                        p(HGNC:MAPK14, pmod(Ph, Tyr, 100))"""
        sub_agent = deepcopy(stmt.enz)
        mc = stmt._get_mod_condition()
        sub_agent.mods.append(mc)
        # FIXME Ignore any bound conditions on the substrate!!!
        # This is because if they are included, a complex node will be returned,
        # which (at least currently) won't incorporate any protein
        # modifications.
        sub_agent.bound_conditions = []
        # FIXME
        self._add_nodes_edges(stmt.enz, sub_agent, pc.DIRECTLY_INCREASES,
                              stmt.evidence)

    def _assemble_transphosphorylation(self, stmt):
        """Example: complex(p(HGNC:EGFR)) =>
                                          p(HGNC:EGFR, pmod(Ph, Tyr, 1173))"""
        # Check our assumptions about the bound condition of the enzyme
        assert len(stmt.enz.bound_conditions) == 1
        assert stmt.enz.bound_conditions[0].is_bound
        # Create a modified protein node for the bound target
        sub_agent = deepcopy(stmt.enz.bound_conditions[0].agent)
        sub_agent.mods.append(stmt._get_mod_condition())
        self._add_nodes_edges(stmt.enz, sub_agent, pc.DIRECTLY_INCREASES,
                              stmt.evidence)

    def _assemble_translocation(self, stmt):
        #cc = hierarchies['cellular_component']
        #nuc_uri = cc.find_entity('nucleus')
        #cyto_uri = cc.find_entity('cytoplasm')
        #cyto_go = cyto_uri.rsplit('/')[-1]
        pass


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


def _get_agent_node(agent):
    if agent.bound_conditions:
        # "Flatten" the bound conditions for the agent at this level
        agent_no_bc = deepcopy(agent)
        agent_no_bc.bound_conditions = []
        members = [agent_no_bc] + [bc.agent for bc in agent.bound_conditions
                                   if bc.is_bound]
        return _get_complex_node(members)
    else:
        return _get_agent_node_no_bcs(agent)


def _get_complex_node(members):
    members_list = []
    for member in members:
        member_data, member_edge = _get_agent_node(member)
        if member_data:
            members_list.append(member_data)
    if members_list:
        complex_node_data = {
                pc.FUNCTION: pc.COMPLEX,
                pc.MEMBERS: members_list}
        return (complex_node_data, None)
    else:
        return (None, None)


def _get_agent_node_no_bcs(agent):
    (abundance_type, db_ns, db_id) = _get_agent_grounding(agent)
    if abundance_type is None:
        logger.warning('Agent %s has no grounding.', agent)
        return (None, None)
    node_data = {pc.FUNCTION: abundance_type,
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
        var = {pc.KIND: pc.HGVS, pc.IDENTIFIER: mut.to_hgvs()}
        variants.append(var)
    if variants:
        node_data[pc.VARIANTS] = variants
    # Also get edge data for the agent
    edge_data = _get_agent_activity(agent)
    return (node_data, edge_data)


def _get_agent_grounding(agent):
    def _get_id(agent, key):
        id = agent.db_refs.get(key)
        if isinstance(id, list):
            id = id[0]
        return id
    hgnc_id = _get_id(agent, 'HGNC')
    uniprot_id = _get_id(agent, 'UP')
    fplx_id = _get_id(agent, 'FPLX')
    ip_id = _get_id(agent, 'IP')
    pfam_id = _get_id(agent, 'PF')
    fa_id = _get_id(agent, 'FA')
    chebi_id = _get_id(agent, 'CHEBI')
    pubchem_id = _get_id(agent, 'PUBCHEM')
    go_id = _get_id(agent, 'GO')
    mesh_id = _get_id(agent, 'MESH')
    if hgnc_id:
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        if not hgnc_name:
            logger.warning('Agent %s with HGNC ID %s has no HGNC name.',
                           agent, hgnc_id)
            return (None, None, None)
        return (pc.PROTEIN, 'HGNC', hgnc_name)
    elif uniprot_id:
        return (pc.PROTEIN, 'UP', uniprot_id)
    elif fplx_id:
        return (pc.PROTEIN, 'FPLX', fplx_id)
    elif pfam_id:
        return (pc.PROTEIN, 'PFAM', pfam_id)
    elif ip_id:
        return (pc.PROTEIN, 'IP', ip_id)
    elif fa_id:
        return (pc.PROTEIN, 'NXPFA', fa_id)
    elif chebi_id:
        if chebi_id.startswith('CHEBI:'):
            chebi_id = chebi_id[len('CHEBI:'):]
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
    text = evidence.text if evidence.text else 'No evidence text.'
    pybel_ev = {pc.EVIDENCE: text}
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

