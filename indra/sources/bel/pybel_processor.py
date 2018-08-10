import re
import logging
from copy import copy
import pybel.constants as pc
from pybel.struct import has_protein_modification
from pybel.canonicalize import edge_to_bel
from indra.statements import *
from indra.sources.bel.belrdf_processor import bel_to_indra, chebi_name_id
from indra.databases import hgnc_client, uniprot_client
from indra.assemblers.pybel_assembler import _pybel_indra_act_map


logger = logging.getLogger('pybel_processor')


_pybel_indra_pmod_map = {
    'Ph': 'phosphorylation',
    'Hy': 'hydroxylation',
    'Sumo': 'sumoylation',
    'Ac': 'acetylation',
    'Glyco': 'glycosylation',
    'ADPRib': 'ribosylation',
    'Ub': 'ubiquitination',
    'Farn': 'farnesylation',
    'Gerger': 'geranylgeranylation',
    'Palm': 'palmitoylation',
    'Myr': 'myristoylation',
    'Me': 'methylation',
}

#: A mapping from the BEL text location annotation to the INDRA ones at
# :py:data:`indra.reach.processor._section_list`
#: see https://arty.scai.fraunhofer.de/artifactory/bel/annotation/text-location/text-location-1.0.0.belanno
_pybel_text_location_map = {
    "Abstract": 'abstract',
    "Results": 'results',
    "Legend": 'figure',
    "Review": None,
    'Introduction': 'introduction',
    'Methods': 'methods',
    'Discussion': 'discussion',
    'Conclusion': 'conclusion'
}


class PybelProcessor(object):
    """Extract INDRA Statements from a PyBEL Graph.

    Currently does not handle non-causal relationships (positiveCorrelation,
    (negativeCorrelation, hasVariant, etc.)

    Parameters
    ----------
    graph : pybel.BELGraph
        PyBEL graph containing the BEL content.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of extracted INDRA Statements representing BEL Statements.
    """
    def __init__(self, graph):
        self.graph = graph
        self.statements = []
        self.unhandled = []

    # FIXME: Handle reactions
    def get_statements(self):
        for u, v, k, d in self.graph.edges_iter(keys=True, data=True):
            u_data = self.graph.node[u]
            v_data = self.graph.node[v]
            # We only interpret causal relations, not correlations
            if d[pc.RELATION] not in pc.CAUSAL_RELATIONS:
                self.unhandled.append((u_data, v_data, k, d))
                continue
            # If the left or right-hand sides involve complex abundances,
            # add them as statements
            for node_ix, node_data in enumerate((u_data, v_data)):
                if node_data[pc.FUNCTION] == pc.COMPLEX:
                    self._get_complex(u_data, v_data, k, d, node_ix)
            subj_activity = _get_activity_condition(d.get(pc.SUBJECT))
            obj_activity = _get_activity_condition(d.get(pc.OBJECT))
            obj_to_loc = _get_translocation_target(d.get(pc.OBJECT))
            # If the object is a translocation, this represents a controlled
            # translocation, which we currently do not represent
            if obj_to_loc:
                self.unhandled.append((u_data, v_data, k, d))
                logger.info("Controlled translocations are currently not "
                            "handled: %s)", edge_to_bel(u_data, v_data, d))
                continue

            v_func = v_data[pc.FUNCTION]
                
            # Modification, e.g.
            #   x(Foo) -> p(Bar, pmod(Ph))
            #   act(x(Foo)) -> p(Bar, pmod(Ph))
            if v_func == pc.PROTEIN and \
               has_protein_modification(self.graph, v):
                if obj_activity:
                    logger.info("Ignoring object activity modifier in "
                                "modification statement: %s, %s, %s, %s",
                                u_data, v_data, k, d)
                else:
                    self._get_modification(u_data, v_data, k, d)
            elif obj_activity:
                # If the agents on the left and right hand sides are the same,
                # then get an active form:
                # ActiveForm
                #   p(Foo, {variants}) ->/-| act(p(Foo))
                # Also Composite active forms:
                #   compositeAbundance(p(Foo, pmod('Ph', 'T')),
                #                       p(Foo, pmod('Ph', 'Y'))) ->/-|
                #                            act(p(Foo))
                if not subj_activity and _proteins_match(u_data, v_data):
                    self._get_active_form(u_data, v_data, k, d)
                # Gef
                #   act(p(Foo)) => gtp(p(Foo))
                # Gap
                #   act(p(Foo)) =| gtp(p(Foo))
                elif subj_activity and _rel_is_direct(d) and \
                     obj_activity.activity_type == 'gtpbound':
                    self._get_gef_gap(u_data, v_data, k, d)
                # Activation/Inhibition
                #   x(Foo) -> act(x(Foo))
                #   act(x(Foo)) -> act(x(Foo))
                # GtpActivation
                #   gtp(p(Foo)) => act(p(Foo))
                else:
                    self._get_regulate_activity(u_data, v_data, k, d)
            # Activations involving biological processes or pathologies
            #   x(Foo) -> bp(Bar)
            elif v_func in (pc.BIOPROCESS, pc.PATHOLOGY):
                self._get_regulate_activity(u_data, v_data, k, d)
            # Regulate amount
            #   x(Foo) -> p(Bar)
            #   x(Foo) -> r(Bar)
            #   act(x(Foo)) -> p(Bar):
            #   x(Foo) -> deg(p(Bar))
            #   act(x(Foo)) ->/-| deg(p(Bar))
            elif v_data[pc.FUNCTION] in (pc.PROTEIN, pc.RNA, pc.ABUNDANCE,
                 pc.COMPLEX, pc.MIRNA) and not obj_activity:
                self._get_regulate_amount(u_data, v_data, k, d)
            # Controlled conversions
            #   x(Foo) -> rxn(reactants(r1,...,rn), products(p1,...pn))
            #   act(x(Foo)) -> rxn(reactants(r1,...,rn), products(p1,...pn))
            # Note that we can't really handle statements where the relation
            # is decreases, as inhibition of a reaction match the semantics
            # of a controlled conversion
            elif v_data[pc.FUNCTION] == pc.REACTION and \
                 d[pc.RELATION] in pc.CAUSAL_INCREASE_RELATIONS:
                self._get_conversion(u_data, v_data, k, d)
            # UNHANDLED
            # rxn(reactants(r1,...,rn), products(p1,...pn))
            # Complex(a,b)
            # p(A, pmod('ph')) -> Complex(A, B)
            # Complex(A-Ph, B) 
            # Complexes
            #   complex(x(Foo), x(Bar), ...)
            else:
                self.unhandled.append((u_data, v_data, k, d))

    def _get_complex(self, u_data, v_data, k, edge_data, node_ix):
        # Get an agent with bound conditions from the Complex
        assert node_ix in (0, 1)
        node_data = [u_data, v_data][node_ix]
        cplx_agent = _get_agent(node_data, None)
        if cplx_agent is None:
            return
        agents = [bc.agent for bc in cplx_agent.bound_conditions]
        cplx_agent.bound_conditions = []
        agents.append(cplx_agent)
        ev = _get_evidence(u_data, v_data, k, edge_data)
        stmt = Complex(agents, evidence=[ev])
        self.statements.append(stmt)

    def _get_regulate_amount(self, u_data, v_data, k, edge_data):
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        obj_agent = _get_agent(v_data, edge_data.get(pc.OBJECT))
        if subj_agent is None or obj_agent is None:
            self.unhandled.append((u_data, v_data, edge_data))
            return
        obj_mod = edge_data.get(pc.OBJECT)
        deg_polarity = (-1 if obj_mod and obj_mod[pc.MODIFIER] == pc.DEGRADATION
                        else 1)
        rel_polarity = (1 if edge_data[pc.RELATION] in
                        pc.CAUSAL_INCREASE_RELATIONS else -1)
        # Set polarity accordingly based on the relation type and whether
        # the object is a degradation node
        if deg_polarity * rel_polarity > 0:
            stmt_class = IncreaseAmount
        else:
            stmt_class = DecreaseAmount
        ev = _get_evidence(u_data, v_data, k, edge_data)
        stmt = stmt_class(subj_agent, obj_agent, evidence=[ev])
        self.statements.append(stmt)

    def _get_modification(self, u_data, v_data, k, edge_data):
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        mods, muts = _get_all_pmods(v_data, edge_data)
        v_data_no_mods = _remove_pmods(v_data)
        obj_agent = _get_agent(v_data_no_mods, edge_data.get(pc.OBJECT))
        if subj_agent is None or obj_agent is None:
            self.unhandled.append((u_data, v_data, k, edge_data))
            return
        for mod in mods:
            modclass = modtype_to_modclass[mod.mod_type]
            ev = _get_evidence(u_data, v_data, k, edge_data)
            stmt = modclass(subj_agent, obj_agent, mod.residue, mod.position,
                            evidence=[ev])
            self.statements.append(stmt)

    def _get_regulate_activity(self, u_data, v_data, k, edge_data):
        # Subject info
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        subj_activity = _get_activity_condition(edge_data.get(pc.SUBJECT))
        subj_function = u_data.get(pc.FUNCTION)
        # Object info
        # Note: Don't pass the object modifier data because we don't want to
        # put an activity on the agent
        obj_agent = _get_agent(v_data, None)
        obj_function = v_data.get(pc.FUNCTION)
        # If it's a bioprocess object, we won't have an activity in the edge
        if obj_function in (pc.BIOPROCESS, pc.PATHOLOGY):
            activity_type = 'activity'
        else:
            obj_activity_condition = \
                            _get_activity_condition(edge_data.get(pc.OBJECT))
            activity_type = obj_activity_condition.activity_type
            assert obj_activity_condition.is_active is True
        # Check for valid subject/object
        if subj_agent is None or obj_agent is None:
            self.unhandled.append((u_data, v_data, edge_data))
            return
        # Check which kind of statement we need to make
        # GtpActivation
        if subj_activity and subj_activity.activity_type == 'gtpbound' and \
           subj_function == pc.PROTEIN and obj_function == pc.PROTEIN and \
           edge_data[pc.RELATION] == pc.DIRECTLY_INCREASES:
            stmt_class = GtpActivation
        elif edge_data[pc.RELATION] in pc.CAUSAL_INCREASE_RELATIONS:
            stmt_class = Activation
        else:
            stmt_class = Inhibition
        ev = _get_evidence(u_data, v_data, k, edge_data)
        stmt = stmt_class(subj_agent, obj_agent, activity_type, evidence=[ev])
        self.statements.append(stmt)

    def _get_active_form(self, u_data, v_data, k, edge_data):
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        # Don't pass the object modifier info because we don't want an activity
        # condition applied to the agent
        obj_agent = _get_agent(v_data)
        if subj_agent is None or obj_agent is None:
            self.unhandled.append((u_data, v_data, edge_data))
            return
        obj_activity_condition = \
                            _get_activity_condition(edge_data.get(pc.OBJECT))
        activity_type = obj_activity_condition.activity_type
        # If the relation is DECREASES, this means that this agent state
        # is inactivating
        is_active = edge_data[pc.RELATION] in pc.CAUSAL_INCREASE_RELATIONS
        ev = _get_evidence(u_data, v_data, k, edge_data)
        stmt = ActiveForm(subj_agent, activity_type, is_active, evidence=[ev])
        self.statements.append(stmt)

    def _get_gef_gap(self, u_data, v_data, k, edge_data):
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        obj_agent = _get_agent(v_data)
        if subj_agent is None or obj_agent is None:
            self.unhandled.append((u_data, v_data, k, edge_data))
            return
        ev = _get_evidence(u_data, v_data, k, edge_data)
        if edge_data[pc.RELATION] in pc.CAUSAL_INCREASE_RELATIONS:
            stmt_class = Gef
        else:
            stmt_class = Gap
        stmt = stmt_class(subj_agent, obj_agent, evidence=[ev])
        self.statements.append(stmt)

    def _get_conversion(self, u_data, v_data, k, edge_data):
        subj_agent = _get_agent(u_data, edge_data.get(pc.SUBJECT))
        # Get the nodes for the reactants and products
        reactant_agents = [_get_agent(r) for r in v_data[pc.REACTANTS]]
        product_agents = [_get_agent(p) for p in v_data[pc.PRODUCTS]]
        if subj_agent is None or \
           any([r is None for r in reactant_agents]) or \
           any([p is None for p in product_agents]):
            self.unhandled.append((u_data, v_data, k, edge_data))
            return
        ev = _get_evidence(u_data, v_data, k, edge_data)
        stmt = Conversion(subj_agent, obj_from=reactant_agents,
                          obj_to=product_agents, evidence=ev)
        self.statements.append(stmt)


def _get_agent(node_data, node_modifier_data=None):
    # FIXME: Handle translocations on the agent for ActiveForms, turn into
    # location conditions
    # Check the node type/function
    node_func = node_data[pc.FUNCTION]
    if node_func not in (pc.PROTEIN, pc.RNA, pc.BIOPROCESS, pc.COMPLEX,
                         pc.PATHOLOGY, pc.ABUNDANCE, pc.MIRNA):
        mod_data = ('No node data' if not node_modifier_data
                                   else node_modifier_data.get(pc.CNAME))
        logger.info("Nodes of type %s not handled: %s",
                    node_func, mod_data)
        return None
    # Skip gene/protein fusions
    if pc.FUSION in node_data:
        logger.info("Gene and protein fusions not handled: %s" % str(node_data))
        return None
    # COMPLEXES ------------
    # First, handle complexes, which will consist recursively of other agents
    if node_func == pc.COMPLEX:
        # First, check for members: if there are no members, we assume this
        # is a named complex
        members = node_data.get(pc.MEMBERS)
        if members is None:
            return None
        # Otherwise, get the "main" agent, to which the other members will be
        # attached as bound conditions
        main_agent = _get_agent(members[0])
        # If we can't get the main agent, return None
        if main_agent is None:
            return None
        bound_conditions = [BoundCondition(_get_agent(m), True)
                            for m in members[1:]]
        # Check the bound_conditions for any None agents
        if any([bc.agent is None for bc in bound_conditions]):
            return None
        main_agent.bound_conditions = bound_conditions
        # Get activity of main agent
        ac = _get_activity_condition(node_modifier_data)
        main_agent.activity = ac
        return main_agent
    # OTHER NODE TYPES -----
    # Get node identifier information
    name = node_data.get(pc.NAME)
    ns = node_data[pc.NAMESPACE]
    ident = node_data.get(pc.IDENTIFIER)
    # No ID present, get identifier using the name, namespace
    db_refs = None
    if not ident:
        assert name, "Node must have a name if lacking an identifier."
        if ns == 'HGNC':
            hgnc_id = hgnc_client.get_hgnc_id(name)
            if not hgnc_id:
                logger.info("Invalid HGNC name: %s (%s)" % (name, node_data))
                return None
            db_refs = {'HGNC': hgnc_id}
            up_id = _get_up_id(hgnc_id)
            if up_id:
                db_refs['UP'] = up_id
        # FIXME: Look up go ID in ontology lookup service
        # FIXME: Look up MESH IDs from name
        # FIXME: For now, just use node name
        elif ns in ('GOBP', 'MESHPP', 'MESHD'):
            db_refs = {}
        # For now, handle MGI/RGD but putting the name into the db_refs so
        # it's clear what namespace the name belongs to
        # FIXME: Full implementation would look up MGI/RGD identifiers from
        # the names, and obtain corresponding Uniprot IDs
        elif ns in ('MGI', 'RGD'):
            db_refs = {ns: name}
        # Map Selventa families to FamPlexes
        elif ns == 'SFAM':
            db_refs = {'SFAM': name}
            indra_name = bel_to_indra.get(name)
            if indra_name is None:
                logger.info('Could not find mapping for BEL/SFAM family: '
                            '%s (%s)' % (name, node_data))
            else:
                db_refs['FPLX'] = indra_name
                name = indra_name
        # Map Entrez genes to HGNC/UP
        elif ns == 'EGID':
            hgnc_id = hgnc_client.get_hgnc_from_entrez(name)
            db_refs = {'EGID': name}
            if hgnc_id is not None:
                db_refs['HGNC'] = hgnc_id
                name = hgnc_client.get_hgnc_name(hgnc_id)
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if up_id:
                    db_refs['UP'] = up_id
                else:
                    logger.info('HGNC entity %s with HGNC ID %s has no '
                                'corresponding Uniprot ID.',
                                name, hgnc_id)
            else:
                logger.info('Could not map EGID%s to HGNC.' % name)
                name = 'E%s' % name
        # CHEBI
        elif ns == 'CHEBI':
            chebi_id = chebi_name_id.get(name)
            if chebi_id:
                db_refs = {'CHEBI': chebi_id}
            else:
                logger.info('CHEBI name %s not found in map.' % name)
        # SDIS, SCHEM: Include the name as the ID for the namespace
        elif ns in ('SDIS', 'SCHEM'):
            db_refs = {ns: name}
        else:
            print("Unhandled namespace: %s: %s (%s)" % (ns, name, node_data))
    # We've already got an identifier, look up other identifiers if necessary
    else:
        # Get the name, overwriting existing name if necessary
        if ns == 'HGNC':
            name = hgnc_client.get_hgnc_name(ident)
            db_refs = {'HGNC': ident}
            up_id = _get_up_id(ident)
            if up_id:
                db_refs['UP'] = up_id
        elif ns == 'UP':
            db_refs = {'UP': ident}
            name = uniprot_client.get_gene_name(ident)
            assert name
            if uniprot_client.is_human(ident):
                hgnc_id = hgnc_client.get_hgnc_id(name)
                if not hgnc_id:
                    logger.info('Uniprot ID linked to invalid human gene '
                                'name %s' % name)
                else:
                    db_refs['HGNC'] = hgnc_id
        elif ns in ('MGI', 'RGD'):
            raise ValueError('Identifiers for MGI and RGD databases are not '
                             'currently handled: %s' % node_data)
        else:
            print("Unhandled namespace with identifier: %s: %s (%s)" %
                  (ns, name, node_data))
    if db_refs is None:
        logger.info('Unable to get identifier information for node: %s'
                     % node_data)
        return None
    # Get modification conditions
    mods, muts = _get_all_pmods(node_data)
    # Get activity condition
    ac = _get_activity_condition(node_modifier_data)
    to_loc = _get_translocation_target(node_modifier_data)
    # Check for unhandled node modifiers, skip if so
    if _has_unhandled_modifiers(node_modifier_data):
        return None
    # Make the agent
    ag = Agent(name, db_refs=db_refs, mods=mods, mutations=muts, activity=ac,
               location=to_loc)
    return ag


def _get_evidence(u_data, v_data, k, edge_data):
    ev_text = edge_data.get(pc.EVIDENCE)
    ev_citation = edge_data.get(pc.CITATION)
    ev_pmid = None
    if ev_citation:
        cit_type = ev_citation[pc.CITATION_TYPE]
        cit_ref = ev_citation[pc.CITATION_REFERENCE]
        if cit_type == pc.CITATION_TYPE_PUBMED:
            ev_pmid = cit_ref
            ev_ref = None
        else:
            ev_pmid = None
            ev_ref = '%s: %s' % (cit_type, cit_ref)
    epistemics = {'direct': _rel_is_direct(edge_data)}
    annotations = edge_data.get(pc.ANNOTATIONS, {})
    annotations['bel'] = edge_to_bel(u_data, v_data, edge_data)
    if ev_ref:  # FIXME what if ev_citation is Falsy?
        annotations['citation_ref'] = ev_ref

    text_location = annotations.pop('TextLocation', None)
    if text_location:
        # Handle dictionary text_location like {'Abstract': True}
        if isinstance(text_location, dict):
            # FIXME: INDRA's section_type entry is meant to contain
            # a single section string like "abstract" but in principle
            # pybel could have a list of entries in the TextLocation dict.
            # Here we just take the first one.
            text_location = list(text_location.keys())[0]
        epistemics['section_type'] = _pybel_text_location_map.get(text_location)

    ev = Evidence(text=ev_text, pmid=ev_pmid, source_api='bel',
                  source_id=k, epistemics=epistemics,
                  annotations=annotations)
    return ev


def _rel_is_direct(d):
    return d[pc.RELATION] in (pc.DIRECTLY_INCREASES, pc.DIRECTLY_DECREASES)


def _get_up_id(hgnc_id):
    up_id = hgnc_client.get_uniprot_id(hgnc_id)
    if not up_id:
        logger.info("No Uniprot ID for HGNC ID %s" % hgnc_id)
    return up_id


def _remove_pmods(node_data):
    node_data_no_pmods = copy(node_data)
    variants = node_data.get(pc.VARIANTS)
    if variants:
        node_data_no_pmods[pc.VARIANTS] = [
            var for var in variants
            if var[pc.KIND] != pc.PMOD
        ]
    return node_data_no_pmods


def _get_all_pmods(node_data, remove_pmods=False):
    mods = []
    muts = []
    variants = node_data.get(pc.VARIANTS)
    if not variants:
        return mods, muts

    for var in variants:
        if var[pc.KIND] == pc.HGVS:
            hgvs_str = var[pc.IDENTIFIER]
            position, res_from, res_to = _parse_mutation(hgvs_str)
            if position is None and res_from is None and res_to is None:
                logger.info("Could not parse HGVS string %s" % hgvs_str)
            else:
                mut_cond = MutCondition(position, res_from, res_to)
                muts.append(mut_cond)
        elif var[pc.KIND] == pc.PMOD:
            var_id_dict = var[pc.IDENTIFIER]
            var_ns = var_id_dict[pc.NAMESPACE]
            if var_ns == pc.BEL_DEFAULT_NAMESPACE:
                var_id = var_id_dict[pc.NAME]
                mod_type = _pybel_indra_pmod_map.get(var_id)
                if mod_type is None:
                    logger.info("Unhandled modification type %s (%s)" %
                                (var_id, node_data))
                    continue
                mc = ModCondition(mod_type, var.get(pc.PMOD_CODE),
                                  var.get(pc.PMOD_POSITION))
                mods.append(mc)
        # FIXME These unhandled mod types should result in throwing out
        # the node (raise, or return None)
        elif var[pc.KIND] == pc.GMOD:
            logger.debug('Unhandled node variant GMOD: %s' % node_data)
        elif var[pc.KIND] == pc.FRAGMENT:
            logger.debug('Unhandled node variant FRAG: %s' % node_data)
        else:
            logger.debug('Unknown node variant type: %s' % node_data)
    return (mods, muts)


def _get_activity_condition(node_modifier_data):
    if node_modifier_data is None or node_modifier_data == {}:
        return None
    if node_modifier_data[pc.MODIFIER] != pc.ACTIVITY:
        return None
    effect = node_modifier_data.get(pc.EFFECT)
    # No specific effect, just return generic activity
    if not effect:
        return ActivityCondition('activity', True)

    activity_ns = effect[pc.NAMESPACE]
    if activity_ns == pc.BEL_DEFAULT_NAMESPACE:
        activity_name = effect[pc.NAME]
        activity_type = _pybel_indra_act_map.get(activity_name)
        # If an activity type in Bel/PyBel that is not implemented in INDRA,
        # return generic activity
        if activity_type is None:
            return ActivityCondition('activity', True)
        return ActivityCondition(activity_type, True)
    # If an unsupported namespace, simply return generic activity
    return ActivityCondition('activity', True)


def _get_translocation_target(node_modifier_data):
    # First check if there is a translocation modifier
    if node_modifier_data is None or node_modifier_data == {}:
        return None
    if node_modifier_data[pc.MODIFIER] != pc.TRANSLOCATION:
        return None
    # Next, make sure there is information on the translocation target
    transloc_data = node_modifier_data[pc.EFFECT]
    to_loc_info = transloc_data.get(pc.TO_LOC)
    if not to_loc_info:
        return None
    to_loc_ns = to_loc_info.get(pc.NAMESPACE)
    to_loc_name = to_loc_info.get(pc.NAME)
    # Only use GO Cellular Component location names
    if to_loc_ns != 'GOCC' or not to_loc_name:
        return None
    try:
        valid_loc = get_valid_location(to_loc_name)
    except InvalidLocationError:
        return None
    return valid_loc


def _has_unhandled_modifiers(node_modifier_data):
    # First check if there is a translocation modifier
    if node_modifier_data is None or node_modifier_data == {}:
        return False
    mod = node_modifier_data.get(pc.MODIFIER)
    if mod is None:
        return False
    if mod in (pc.CELL_SECRETION, pc.CELL_SURFACE_EXPRESSION):
        logger.info("Unhandled node modifier data: %s" % node_modifier_data)
        return True


def _proteins_match(u_data, v_data):
    return (
        u_data[pc.FUNCTION] == pc.PROTEIN and
        v_data[pc.FUNCTION] == pc.PROTEIN and
        pc.NAMESPACE in u_data and pc.NAMESPACE in v_data and
        pc.NAME in u_data and pc.NAME in v_data and
        u_data[pc.NAMESPACE] == v_data[pc.NAMESPACE] and
        u_data[pc.NAME] == v_data[pc.NAME]
    )


_hgvs_protein_mutation = re.compile('^p.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})')


def _parse_mutation(s):
    m = _hgvs_protein_mutation.match(s)
    if not m:
        return (None, None, None)
    from_aa, position, to_aa = m.groups()
    return position, from_aa, to_aa

