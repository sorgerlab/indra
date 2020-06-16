import re
import copy
import logging
import itertools
import collections
import pybiopax.biopax as bp
from functools import lru_cache

from pybiopax import model_to_owl_file

from indra.databases import hgnc_client, uniprot_client, chebi_client
from indra.statements import *
from indra.util import decode_obj, flatten

logger = logging.getLogger(__name__)

# TODO:
# - Extract cellularLocation from each PhysicalEntity
# - Extract and represent FragmentFeatures
# - Extract and represent BindingFeatures
# - direct is not always appropriate, e.g., for CTD as a source
# - Look at participantStoichiometry within BiochemicalReaction
# - Check whether to use Control or only Catalysis (Control might not
#   be direct)
# - Implement extracting modifications with Complex enzyme
# - Implement extracting modifications with Complex substrate


class BiopaxProcessor(object):
    """The BiopaxProcessor extracts INDRA Statements from a BioPAX model.

    The BiopaxProcessor uses pattern searches in a BioPAX OWL model to
    extract mechanisms from which it constructs INDRA Statements.

    Parameters
    ----------
    model : org.biopax.paxtools.model.Model
        A BioPAX model object (java object)

    Attributes
    ----------
    model : org.biopax.paxtools.model.Model
        A BioPAX model object (java object) which is queried using Paxtools
        to extract INDRA Statements
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted from the model.
    """
    def __init__(self, model):
        self.model = model
        self.statements = []
        self._mod_conditions = {}
        self._activity_conditions = {}
        self._agents = {}

    def print_statements(self):
        """Print all INDRA Statements collected by the processors."""
        for i, stmt in enumerate(self.statements):
            print("%s: %s" % (i, stmt))

    def save_model(self, file_name):
        """Save the BioPAX model object in an OWL file.

        Parameters
        ----------
        file_name : str
            The name of the OWL file to save the model in.
        """
        model_to_owl_file(self.model, file_name)

    def eliminate_exact_duplicates(self):
        """Eliminate Statements that were extracted multiple times.

        Due to the way the patterns are implemented, they can sometimes yield
        the same Statement information multiple times, in which case,
        we end up with redundant Statements that aren't from independent
        underlying entries. To avoid this, here, we filter out such
        duplicates.
        """
        # Here we use the deep hash of each Statement, and by making a dict,
        # we effectively keep only one Statement with a given deep hash
        self.statements = list({stmt.get_hash(shallow=False, refresh=True): stmt
                                for stmt in self.statements}.values())

    @staticmethod
    def find_matching_left_right(conversion):
        matches = []

        left_simple = conversion.left[:]
        for pe in left_simple:
            if _is_complex(pe):
                left_simple += pe.component
        right_simple = conversion.right[:]
        for pe in right_simple:
            if _is_complex(pe):
                right_simple += pe.component

        for inp, outp in itertools.product(left_simple, right_simple):
            if _is_simple_physical_entity(inp) and \
                    _is_simple_physical_entity(outp):
                if inp.entity_reference == outp.entity_reference:
                    matches.append((inp, outp))
        return matches

    def feature_delta(self, from_pe: bp.PhysicalEntity,
                      to_pe: bp.PhysicalEntity):
        # First deal with activity changes
        from_acts = {self._activity_conditions[f.uid] for f in from_pe.feature
                     if f.uid in self._activity_conditions}
        to_acts = {self._activity_conditions[f.uid] for f in to_pe.feature
                   if f.uid in self._activity_conditions}
        assert len(from_acts) <= 1
        assert len(to_acts) <= 1
        activity_change = None
        if from_acts == {'active'}:
            if not to_acts or to_acts == {'inactive'}:
                activity_change = 'inactive'
        elif not from_acts or from_acts == {'inactive'}:
            if to_acts == {'active'}:
                activity_change = 'active'

        # Now look at modification changes
        from_mods = {self._mod_conditions[f.uid] for f in from_pe.feature
                     if f.uid in self._mod_conditions}
        from_mods = {mc.matches_key(): mc for mc in from_mods}
        to_mods = {self._mod_conditions[f.uid] for f in to_pe.feature
                   if f.uid in self._mod_conditions}
        to_mods = {mc.matches_key(): mc for mc in to_mods}

        gained_mods = {to_mods[k] for k in
                       set(to_mods.keys()) - set(from_mods.keys())}
        lost_mods = {from_mods[k] for k in
                     set(from_mods.keys()) - set(to_mods.keys())}
        return gained_mods, lost_mods, activity_change

    def extract_features(self):
        for feature in self.model.get_objects_by_type(bp.EntityFeature):
            # TODO: handle BindingFeatures and FragmentFeatures here
            if not isinstance(feature, bp.ModificationFeature):
                continue
            mf_type = get_modification_type(feature)
            if mf_type == 'modification':
                mod = self.mod_condition_from_mod_feature(feature)
                if mod:
                    self._mod_conditions[feature.uid] = mod
            elif mf_type == 'activity':
                self._activity_conditions[feature.uid] = 'active'
            elif mf_type == 'inactivity':
                self._activity_conditions[feature.uid] = 'inactive'

    def get_complexes(self):
        """Extract INDRA Complex Statements from the BioPAX model.

        This method searches for org.biopax.paxtools.model.level3.Complex
        objects which represent molecular complexes. It doesn't reuse
        BioPAX Pattern's org.biopax.paxtools.pattern.PatternBox.inComplexWith
        query since that retrieves pairs of complex members rather than
        the full complex.
        """
        for bpe in self.model.get_objects_by_type(bp.Complex):
            ev = self._get_evidence(bpe)
            members = self._get_complex_members(bpe)
            if members is not None:
                if len(members) > 10:
                    logger.debug('Skipping complex with more than 10 members.')
                    continue
                complexes = _get_combinations(members)
                for c in complexes:
                    self.statements.append(Complex(c, ev))

    def _conversion_state_iter(self):
        for control in self.model.get_objects_by_type(bp.Control):
            conversion = control.controlled
            # Sometimes there is nothing being controlled, we skip
            # those cases. We also want to skip things like Modulation
            # that don't convert anything.
            if not isinstance(conversion, bp.Conversion):
                continue
            control_agents = flatten([self._get_primary_controller(c) for c in
                                      control.controller])
            control_agents = [c for c in control_agents if c is not None]
            if not control_agents:
                continue
            ev = self._get_evidence(control)
            for inp, outp in self.find_matching_left_right(conversion):
                inp_agents = self._get_agents_from_entity(inp)
                gained_mods, lost_mods, activity_change = \
                    self.feature_delta(inp, outp)
                yield control_agents, inp_agents, gained_mods, \
                    lost_mods, activity_change, ev

    def get_modifications(self):
        """Extract INDRA Modification Statements from the BioPAX model."""
        for control_agents, inp_agents, gained_mods, lost_mods, \
                activity_change, ev in self._conversion_state_iter():
            for mods, is_gain in ((gained_mods, True), (lost_mods, False)):
                for mod in mods:
                    stmt_class = modtype_to_modclass[mod.mod_type]
                    if not is_gain:
                        stmt_class = modclass_to_inverse[stmt_class]
                    for enz, sub in \
                            itertools.product(_listify(control_agents),
                                              _listify(inp_agents)):
                        assert isinstance(enz, Agent)
                        assert isinstance(sub, Agent)
                        stmt = stmt_class(enz, sub, mod.residue,
                                          mod.position, evidence=ev)
                        stmt = _remove_redundant_mods(stmt)
                        self.statements.append(stmt)

    def get_regulate_activities(self):
        """Get Activation/Inhibition INDRA Statements from the BioPAX model."""
        for control_agents, inp_agents, gained_mods, lost_mods, \
                activity_change, ev in self._conversion_state_iter():
            # We don't want to have gained or lost modification features
            if not activity_change or (gained_mods or lost_mods):
                continue
            stmt_class = Activation if activity_change == 'active' \
                else Inhibition
            for subj, obj in \
                    itertools.product(_listify(control_agents),
                                      _listify(inp_agents)):
                assert isinstance(subj, Agent)
                assert isinstance(obj, Agent)
                stmt = stmt_class(subj, obj, 'activity',
                                  evidence=ev)
                self.statements.append(stmt)

    def get_activity_modification(self):
        """Extract INDRA ActiveForm statements from the BioPAX model."""
        for control_agents, inp_agents, gained_mods, lost_mods, \
                activity_change, ev in self._conversion_state_iter():
            # We have to have both a modification change and an activity
            # change
            if not (gained_mods or lost_mods) or not activity_change:
                continue
            is_active = (activity_change == 'active')
            for agent in _listify(inp_agents):
                # NOTE: with the ActiveForm representation we cannot
                # separate static_mods and gained_mods. We assume here
                # that the static_mods are inconsequential and therefore
                # are not mentioned as an Agent condition, following
                # don't care don't write semantics. Therefore only the
                # gained_mods are listed in the ActiveForm as Agent
                # conditions.
                assert isinstance(agent, Agent)
                if gained_mods:
                    ag = copy.deepcopy(agent)
                    ag.mods = gained_mods
                    stmt = ActiveForm(ag, 'activity', is_active, evidence=ev)
                    self.statements.append(stmt)
                if lost_mods:
                    ag = copy.deepcopy(agent)
                    ag.mods = lost_mods
                    stmt = ActiveForm(ag, 'activity', not is_active, evidence=ev)
                    self.statements.append(stmt)

    def get_regulate_amounts(self):
        """Extract INDRA RegulateAmount Statements from the BioPAX model."""
        for control in self.model.get_objects_by_type(bp.Control):
            if not isinstance(control.controlled, bp.TemplateReaction):
                continue
            temp_react = control.controlled
            ev = self._get_evidence(control)
            stmt_type = IncreaseAmount if control.control_type == 'ACTIVATION' \
                else DecreaseAmount
            control_agents = flatten([self._get_primary_controller(c) for c in
                                      control.controller])
            if not control_agents:
                continue
            ev = self._get_evidence(control)
            for product in temp_react.product:
                product_agents = self._get_agents_from_entity(product)
                for subj, obj in itertools.product(_listify(control_agents),
                                                   _listify(product_agents)):
                    stmt = stmt_type(subj, obj, evidence=ev)
                    self.statements.append(stmt)

    def get_conversions(self):
        """Extract Conversion INDRA Statements from the BioPAX model.

        This method uses a custom BioPAX Pattern
        (one that is not implemented PatternBox) to query for
        BiochemicalReactions whose left and right hand sides are collections
        of SmallMolecules. This pattern thereby extracts metabolic
        conversions as well as signaling processes via small molecules
        (e.g. lipid phosphorylation or cleavage).
        """
        for control in self.model.get_objects_by_type(bp.Control):
            conversion = control.controlled
            # Sometimes there is nothing being controlled, we skip
            # those cases. We also want to skip things like Modulation
            # that don't convert anything.
            if not isinstance(conversion, bp.Conversion):
                continue

            # We only extract conversions for small molecules
            if not all(_is_small_molecule(pe)
                       for pe in (conversion.left + conversion.right)):
                continue

            # Get the control agents
            control_agents = flatten([self._get_primary_controller(c) for c in
                                      control.controller])
            control_agents = [c for c in control_agents if c is not None]
            if not control_agents:
                continue

            # Assemble from and to object lists
            obj_from = []
            obj_to = []
            for participants, obj_list in ((conversion.left, obj_from),
                                           (conversion.right, obj_to)):
                for participant in participants:
                    agent = self._get_agents_from_entity(participant)
                    if isinstance(agent, list):
                        obj_list += agent
                    else:
                        obj_list.append(agent)

            # Make statements
            ev = self._get_evidence(control)
            for subj in control_agents:
                st = Conversion(subj, obj_from, obj_to, evidence=ev)
                self.statements.append(st)
        return

    def find_gdp_gtp_complex(self, cplxes):
        for cplx in cplxes:
            members = self._get_complex_members(cplx)
            if not members:
                continue
            gdp_gtp_idx = None
            ras_agent = None
            for idx, member in enumerate(members):
                if isinstance(member, Agent) \
                        and member.name in {'GDP', 'GTP'}:
                    gdp_gtp_idx = idx
                    break
            for idx, member in enumerate(members):
                if isinstance(member, Agent) \
                        and 'HGNC' in member.db_refs:
                    ras_agent = member
            if gdp_gtp_idx is None or ras_agent is None:
                continue
            return ras_agent, members[gdp_gtp_idx].name
        return None, None

    def get_gap_gef(self):
        for control in self.model.get_objects_by_type(bp.Control):
            conversion = control.controlled
            # Sometimes there is nothing being controlled, we skip
            # those cases. We also want to skip things like Modulation
            # that don't convert anything.
            if not isinstance(conversion, bp.Conversion):
                continue
            ev = self._get_evidence(control)
            left_complexes = [bpe for bpe in conversion.left
                              if _is_complex(bpe)]
            right_complexes = [bpe for bpe in conversion.right
                               if _is_complex(bpe)]
            left_ras, left_gtp_gdp = \
                self.find_gdp_gtp_complex(left_complexes)
            right_ras, right_gtp_gdp = \
                self.find_gdp_gtp_complex(right_complexes)
            if left_gtp_gdp == 'GDP' and right_gtp_gdp == 'GTP':
                stmt_type = Gef
            elif left_gtp_gdp == 'GTP' and right_gtp_gdp == 'GDP':
                stmt_type = Gap
            else:
                continue

            # Get the control agents
            control_agents = flatten([self._get_primary_controller(c) for c in
                                      control.controller])
            control_agents = [c for c in control_agents if c is not None]
            if not control_agents:
                continue

            ev = self._get_evidence(control)
            for gap_gef, ras in itertools.product(_listify(control_agents),
                                                  _listify(left_ras)):
                st = stmt_type(gap_gef, ras, evidence=ev)
                self.statements.append(st)

    def _get_complex_members(self, cplx: bp.Complex):
        # Get the members of a complex. This is returned as a list
        # of lists since complexes can contain other complexes. The
        # list of lists solution allows us to preserve this.
        member_pes = cplx.component

        # Make a dict of member URIs and their
        # corresponding stoichiometries
        member_stos = {cs.physical_entity.uid: cs.stoichiometric_coefficient
                       for cs in cplx.component_stoichiometry}

        # Some complexes do not have any members explicitly listed
        if not member_pes:
            member_pes = cplx.member_physical_entity
            if not member_pes:
                logger.debug('Complex "%s" has no members.' %
                             cplx.display_name)
                return None
        members = []
        for m in member_pes:
            if _is_complex(m):
                ms = self._get_complex_members(m)
                if ms is None:
                    return None
                members.extend(ms)
            else:
                ma = self._get_agents_from_entity(m)
                try:
                    sto = member_stos[m.uid]
                    sto_int = int(float(sto))  # This is needed for e.g., '1.0'
                except KeyError:
                    # No stoichiometry information - assume it is 1
                    sto_int = 1
                for i in range(sto_int):
                    members.append(ma)
        return members

    @staticmethod
    def _get_entity_mods(bpe):
        """Get all the modifications of an entity in INDRA format"""
        if _is_entity(bpe):
            features = bpe.feature
        else:
            features = bpe.entity_feature
        mods = []
        for feature in features:
            if not _is_modification(feature):
                continue
            mc = BiopaxProcessor.mod_condition_from_mod_feature(feature)
            if mc is not None:
                mods.append(mc)
        return mods

    def _get_primary_controller(self, controller_pe):
        if not isinstance(controller_pe, bp.PhysicalEntity):
            return None

        # If it's not a complex, just return the corresponding agent
        if not _is_complex(controller_pe):
            enzs = self._get_agents_from_entity(controller_pe)
            return enzs

        # Identifying the "real" enzyme in a complex may not always be
        # possible.
        # One heuristic here could be to find the member which is
        # active and if it is the only active member then
        # set this as the enzyme to which all other members of the
        # complex are bound.
        # Get complex members
        members = self._get_complex_members(controller_pe)
        if members is None:
            return None
        # Separate out protein and non-protein members
        protein_members = []
        non_protein_members = []
        for m in members:
            if isinstance(m, Agent):
                if m.db_refs.get('UP') or \
                        m.db_refs.get('HGNC'):
                    protein_members.append(m)
                else:
                    non_protein_members.append(m)
            else:
                all_protein = True
                for subm in m:
                    if not (subm.db_refs.get('UP') or \
                            subm.db_refs.get('HGNC')):
                        all_protein = False
                        break
                if all_protein:
                    protein_members.append(m)
                else:
                    non_protein_members.append(m)
        # If there is only one protein member, we can assume that
        # it is the enzyme, and everything else is just bound
        # to it.
        if len(protein_members) == 1:
            enzs = protein_members[0]
            # Iterate over non-protein members
            for bound in non_protein_members:
                if isinstance(bound, Agent):
                    bc = BoundCondition(bound, True)
                    if isinstance(enzs, Agent):
                        enzs.bound_conditions.append(bc)
                    else:
                        for enz in enzs:
                            enz.bound_conditions.append(bc)
                else:
                    msg = 'Cannot handle complex enzymes with ' + \
                            'aggregate non-protein binding partners.'
                    logger.debug(msg)
                    continue
            return enzs
        else:
            msg = 'Cannot handle complex enzymes with ' + \
                    'multiple protein members.'
            logger.debug(msg)
            return None

    def _get_agent_from_entity(self, bpe):
        try:
            return copy.deepcopy(self._agents[bpe.uid])
        except KeyError:
            pass
        name = BiopaxProcessor._get_element_name(bpe)
        db_refs = BiopaxProcessor._get_db_refs(bpe)
        if _is_protein(bpe):
            mcs = BiopaxProcessor._get_entity_mods(bpe)
        else:
            mcs = []
        agent = Agent(name, db_refs=db_refs, mods=mcs)
        return agent

    def _get_agents_from_entity(self, bpe: bp.PhysicalEntity,
                                expand_pe=True, expand_er=True):
        # If the entity has members (like a protein family),
        # we iterate over them
        if expand_pe:
            members = bpe.member_physical_entity
            if members:
                agents = []
                for m in members:
                    member_agents = self._get_agents_from_entity(m)
                    if isinstance(member_agents, Agent):
                        agents.append(member_agents)
                    else:
                        agents.extend(member_agents)
                return agents

        # If the entity has a reference which has members, we iterate
        # over them.
        if expand_er:
            er = BiopaxProcessor._get_entref(bpe)
            if er is not None:
                members = er.member_entity_reference
                if members:
                    agents = []
                    for m in members:
                        agent = self._get_agent_from_entity(m)
                        # For entity references, we remove context
                        agent.mods = []
                        agents.append(agent)
                    return agents
        # If it is a single entity, we get its name and database
        # references
        agent = self._get_agent_from_entity(bpe)
        return agent

    @staticmethod
    def mod_condition_from_mod_feature(mf: bp.ModificationFeature):
        """Extract the type of modification and the position from
        a ModificationFeature object in the INDRA format."""
        # ModificationFeature / SequenceModificationVocabulary
        mf_type = mf.modification_type
        if mf_type is None:
            return None
        for t in mf_type.term:
            if t.startswith('MOD_RES '):
                t = t[8:]
            mf_type_indra = _mftype_dict.get(t)
            if mf_type_indra is not None:
                break
        else:
            logger.debug('Skipping modification with unknown terms: %s' %
                         ', '.join(mf_type.term))
            return None

        mod_type, residue = mf_type_indra

        if mf.feature_location is not None:
            # If it is not a SequenceSite we can't handle it
            if not isinstance(mf.feature_location, bp.SequenceSite):
                mod_pos = None
            else:
                mf_pos_status = mf.feature_location.position_status
                if mf_pos_status is None:
                    mod_pos = None
                elif mf_pos_status and mf_pos_status != 'EQUAL':
                    logger.debug('Modification site position is %s' %
                                 mf_pos_status)
                    mod_pos = None
                else:
                    mod_pos = str(mf.feature_location.sequence_position)
        else:
            mod_pos = None
        mc = ModCondition(mod_type, residue, mod_pos, True)
        return mc

    @staticmethod
    def _get_evidence(bpe: bp.PhysicalEntity):
        citations = BiopaxProcessor._get_citations(bpe)
        if not citations:
            citations = [None]
        epi = {'direct': True}
        annotations = {}
        if bpe.data_source:
            if len(bpe.data_source) > 1:
                logger.warning('More than one data source for %s' % bpe.uid)
            db_name = bpe.data_source[0].display_name
            if db_name:
                annotations['source_sub_id'] = db_name.lower()
        ev = [Evidence(source_api='biopax', pmid=cit,
                       source_id=bpe.uid, epistemics=epi,
                       annotations=annotations)
              for cit in citations]
        return ev

    @staticmethod
    def _get_citations(bpe: bp.PhysicalEntity):
        refs = []
        for xr in bpe.xref:
            db_name = xr.db
            if db_name is not None and db_name.upper() == 'PUBMED':
                refs.append(xr.id)
        # TODO: handle non-pubmed evidence
        # TODO: do we need to look at bpe.getEvidence()
        return refs

    @staticmethod
    def _get_db_refs(bpe: bp.PhysicalEntity):
        db_refs = {}
        if _is_protein(bpe) or _is_rna(bpe):
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            # Handle missing HGNC/UP ids
            if hgnc_id and not uniprot_id:
                uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
            elif uniprot_id and not hgnc_id:
                uniprot_id_lookup = uniprot_id.split('-')[0]
                hgnc_id = uniprot_client.get_hgnc_id(uniprot_id_lookup)
            # If we have both an HGNC ID and a Uniprot ID, override the
            # Uniprot ID with the one associated with the HGNC ID
            elif uniprot_id and hgnc_id:
                hgnc_up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if hgnc_up_id and ',' in hgnc_up_id:
                    up_ids = hgnc_up_id.split(', ')
                    if uniprot_id not in up_ids:
                        logger.info('Uniprot ID %s does not match %s obtained '
                                    'from HGNC ID %s' %
                                    (uniprot_id, hgnc_up_id, hgnc_id))
                elif hgnc_up_id != uniprot_id:
                    logger.info('Uniprot ID %s does not match %s obtained '
                                'from HGNC ID %s' %
                                (uniprot_id, hgnc_up_id, hgnc_id))
            if hgnc_id is not None:
                db_refs['HGNC'] = hgnc_id
            if uniprot_id is not None:
                db_refs['UP'] = uniprot_id
            if not hgnc_id and not uniprot_id:
                rna_groundings = BiopaxProcessor._get_rna_grounding(bpe)
                db_refs.update(rna_groundings)
        elif _is_small_molecule(bpe):
            chebi_id = BiopaxProcessor._get_chebi_id(bpe)
            if chebi_id is not None:
                db_refs['CHEBI'] = chebi_id
            else:
                chemical_groundings = \
                    BiopaxProcessor._get_chemical_grounding(bpe)
                db_refs.update(chemical_groundings)
        else:
            chebi_id = BiopaxProcessor._get_chebi_id(bpe)
            if chebi_id is not None:
                db_refs['CHEBI'] = chebi_id
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            if hgnc_id is not None:
                db_refs['HGNC'] = hgnc_id
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            if uniprot_id is not None:
                db_refs['UP'] = uniprot_id
        return db_refs

    @staticmethod
    @lru_cache(maxsize=1000)
    def _get_element_name(bpe: bp.PhysicalEntity):
        def get_name(bpe):
            # FIXME Deal with case when HGNC entry is not name
            # Deal with case when multiple Uniprot IDs marked as
            # primary
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            if hgnc_id is not None:
                name = hgnc_client.get_hgnc_name(hgnc_id)
                if name is None:
                    name = bpe.display_name
            elif uniprot_id is not None:
                try:
                    name = uniprot_client.get_gene_name(uniprot_id)
                except Exception:
                    name = None
                if name is None:
                    name = bpe.display_name
            else:
                name = bpe.display_name
            return name

        if _is_protein(bpe) or _is_rna(bpe):
            name = get_name(bpe)
        elif _is_small_molecule(bpe):
            name = bpe.display_name
        elif _is_physical_entity(bpe):
            name = bpe.display_name
        else:
            logger.debug('Unhandled entity type %s' %
                         bpe.__class__.__name__)
            name = bpe.display_name

        return name

    @staticmethod
    def _get_uniprot_id(bpe: bp.PhysicalEntity):
        # There is often more than one UniProt ID reported.
        # This usually corresponds to the primary accession ID and one or more
        # secondary accession IDs (these IDs are from deprecated entries that
        # have been merged into the primary).
        def map_to_up_primary(ids):
            primary_ids = set()
            for up_id in ids:
                if not uniprot_client.is_secondary(up_id):
                    primary_ids.add(up_id)
                    continue
                primary_id = uniprot_client.get_primary_id(up_id)
                primary_ids.add(primary_id)
            primary_ids = sorted(list(primary_ids))
            # If there are no primary IDs, we return None
            if not primary_ids:
                return None
            # Try to get primary IDs if there are 
            # If there is more than one primary ID then we return the first one
            elif len(primary_ids) > 1:
                human_upids = [id for id in primary_ids
                               if uniprot_client.is_human(id)]
                if not human_upids:
                    logger.info('More than one primary id but none human, '
                                'choosing the first: %s' %
                                ','.join(primary_ids))
                    primary_id = primary_ids[0]
                elif len(human_upids) > 1:
                    logger.info('More than one human primary id, choosing '
                                'the first: %s' % ','.join(human_upids))
                    primary_id = human_upids[0]
                # Only one, so use it
                else:
                    primary_id = human_upids[0]
            # One primary ID, so use it
            else:
                primary_id = primary_ids[0]
            # Make sure it's unicode
            return str(primary_id)

        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        # First try to match the URI itself to see if it is a UniProt
        # reference.
        m = re.match('http://identifiers.org/uniprot/([A-Z0-9]+)',
                     bp_entref.uid)
        if m:
            uniprot_id = m.groups()[0]
            primary_id = map_to_up_primary([uniprot_id])
            return primary_id
        # If the URI is not a UniProt reference then we look through xrefs
        xrefs = bp_entref.xref
        uniprot_refs = [x for x in xrefs if
                        (x.db is not None and
                         x.db.lower() in ('uniprot knowledgebase',
                                          'uniprotkb',
                                          'uniprot isoform'))]
        if not uniprot_refs:
            return None
        uniprot_ids = [r.id for r in uniprot_refs]
        primary_id = map_to_up_primary(uniprot_ids)
        return primary_id

    @staticmethod
    def _get_chebi_id(bpe: bp.PhysicalEntity):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.xref
        chebi_ids = []
        for xr in xrefs:
            dbname = xr.db
            dbid = xr.id
            if dbname is None:
                continue
            dbname = dbname.upper()
            if dbname == 'CHEBI':
                chebi_ids.append(dbid.replace('CHEBI:', ''))
            elif dbname == 'CAS':
                chebi_mapped = chebi_client.get_chebi_id_from_cas(dbid)
                if chebi_mapped is not None:
                    chebi_ids.append(chebi_mapped)
                else:
                    logger.info('Unknown CAS id: %s (%s)' %
                                 (dbid, bpe.display_name))
        if not chebi_ids:
            return None
        elif len(chebi_ids) == 1:
            return chebi_ids[0]
        else:
            name = BiopaxProcessor._get_element_name(bpe)
            specific_chebi_id = get_specific_chebi_id(frozenset(chebi_ids),
                                                      name)
            return specific_chebi_id

    @staticmethod
    def _get_rna_grounding(bpe: bp.PhysicalEntity):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return {}
        xrefs = bp_entref.xref
        rna_grounding = {}
        for xr in xrefs:
            dbname = xr.db
            dbid = xr.id
            if dbname is None:
                continue
            dbname = dbname.upper()
            if dbname in ('MIRBASE SEQUENCE', 'MIRBASE'):
                rna_grounding['MIRBASE'] = dbid
            elif dbname == 'MIRBASE MATURE SEQUENCE':
                rna_grounding['MIRBASEM'] = dbid
            elif dbname == 'NCBI GENE':
                rna_grounding['NCBI'] = dbid
            elif dbname == 'ENSEMBL':
                rna_grounding['ENSEMBL'] = dbid
        return rna_grounding

    @staticmethod
    def _get_chemical_grounding(bpe: bp.PhysicalEntity):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return {}
        xrefs = bp_entref.xref
        chemical_grounding = {}
        for xr in xrefs:
            dbname = xr.db
            dbid = xr.id
            if dbname is None:
                continue
            dbname = dbname.upper()
            if dbname == 'PUBCHEM-COMPOUND':
                chemical_grounding['PUBCHEM'] = '%s' % dbid
            elif dbname == 'MESH':
                chemical_grounding['MESH'] = dbid
            elif dbname == 'DRUGBANK':
                chemical_grounding['DRUGBANK'] = dbid
            elif dbname == 'HMDB':
                chemical_grounding['HMDB'] = dbid
        return chemical_grounding

    @staticmethod
    def _get_hgnc_id(bpe: bp.PhysicalEntity):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.xref
        # Check for HGNC IDs
        hgnc_ids = [x.id for x in xrefs if
                    (x.db is not None and x.db.lower() == 'hgnc')]
        hgnc_id = None
        for hgnc_id in hgnc_ids:
            m = re.match('([0-9]+)', hgnc_id)
            if m:
                hgnc_id = str(m.groups()[0])
            else:
                m = re.match('hgnc:([0-9]+)', hgnc_id.lower())
                if m:
                    hgnc_id = str(m.groups()[0])
        # If there is no HGNC ID, check for an HGNC symbol and convert back
        # to HGNC
        if not hgnc_id:
            hgnc_syms = [x.id for x in xrefs
                         if (x.db is not None and
                             x.db.lower() == 'hgnc symbol')]
            # If no symbol and no ID, return None
            if not hgnc_syms:
                return None
            # On the off chance that there is more than one symbol, issue
            # a log message and choose the first
            else:
                if len(hgnc_syms) > 1:
                    logger.info('No HGNC ID, and more than one HGNC symbol '
                                'found, using 1st: %s' % str(hgnc_syms))
                hgnc_sym = hgnc_syms[0]
                hgnc_id = hgnc_client.get_hgnc_id(hgnc_sym)
        return hgnc_id

    @staticmethod
    def _get_entref(bpe: bp.PhysicalEntity):
        """Returns the entity reference of an entity if it exists or
        return the entity reference that was passed in as argument."""
        if not _is_reference(bpe):
            try:
                er = bpe.entity_reference
            except AttributeError:
                return None
            return er
        else:
            return bpe

    def get_coverage(self):
        uids = set()
        objs = self.model.objects
        for uid, obj in objs.items():
            if isinstance(obj, (bp.Catalysis, bp.TemplateReactionRegulation)):
                uids.add(uid)
        stmt_uids = set()
        for stmt in self.statements:
            for ev in stmt.evidence:
                stmt_uids.add(ev.source_id)

        uids_not_covered = uids - stmt_uids
        print('Total in model: %d' % len(uids))
        print('Total covered: %d' % len(uids & stmt_uids))
        print('%.2f%% coverage' % (100.0*len(uids & stmt_uids)/len(uids)))
        return len(uids), len(uids & stmt_uids)


_mftype_dict = {
    'phosres': ('phosphorylation', None),
    'phosphorylation': ('phosphorylation', None),
    'phosphorylated residue': ('phosphorylation', None),
    'phosphorylated': ('phosphorylation', None),
    'O-phospho-L-serine': ('phosphorylation', 'S'),
    'O-phosphopantetheine-L-serine': ('phosphorylation', 'S'),
    'opser': ('phosphorylation', 'S'),
    'O-phospho-L-threonine': ('phosphorylation', 'T'),
    'opthr': ('phosphorylation', 'T'),
    'O-phospho-L-tyrosine': ('phosphorylation', 'Y'),
    'O4\'-phospho-L-tyrosine': ('phosphorylation', 'Y'),
    'optyr': ('phosphorylation', 'Y'),
    'ubiquitinated lysine': ('ubiquitination', 'K'),
    'N4-glycosyl-L-asparagine': ('glycosylation', 'N'),
    'n4glycoasn': ('glycosylation', 'N'),
    'O-glycosyl-L-threonine': ('glycosylation', 'T'),
    'S-palmitoyl-L-cysteine': ('palmitoylation', 'C'),
    'N6-acetyllysine': ('acetylation', 'K'),
    'N6-acetyl-L-lysine' : ('acetylation', 'K'),
    'n6aclys': ('acetylation', 'K'),
    'naclys': ('acetylation', 'K'),
    'N-acetylated L-lysine': ('acetylation', 'K'),
    'N-acetylglycine': ('acetylation', 'G'),
    'N-acetylmethionine': ('acetylation', 'M'),
    'Hydroxyproline': ('hydroxylation', 'P'),
    'hydroxylated proline': ('hydroxylation', 'P'),
    '3-hydroxyproline': ('hydroxylation', 'P'),
    '4-hydroxyproline': ('hydroxylation', 'P'),
    '5-hydroxylysine': ('hydroxylation', 'K'),
    'N-myristoylglycine': ('myristoylation', 'G'),
    'N-myristoyl-glycine': ('myristoylation', 'G'),
    'sumoylated lysine': ('sumoylation', 'K'),
    'mearg': ('methylation', 'R'),
    'methylated L-arginine': ('methylation', 'R'),
    'methylated arginine': ('methylation', 'R'),
    'melys' : ('methylation', 'K'),
    'methylated lysine' : ('methylation', 'K'),
    'methylated L-lysine' : ('methylation', 'K'),
    'ubiquitination': ('ubiquitination', None),
    'ubiquitinylated lysine': ('ubiquitination', 'K'),
    'ubiquitination signature tetrapeptidyl lysine': ('ubiquitination', 'K'),
    'Phosphoserine': ('phosphorylation', 'S'),
    'Phosphothreonine': ('phosphorylation', 'T'),
    'Phosphotyrosine': ('phosphorylation', 'Y'),
    'N-acetylalanine': ('acetylation', 'A'),
    'N-acetylserine': ('acetylation', 'S'),
    'N-acetylthreonine': ('acetylation', 'T'),
    'N-acetylvaline': ('acetylation', 'V'),
    'Omega-N-methylarginine': ('methylation', 'R'),
    'N6-methyllysine': ('methylation', 'K'),
    'Dimethylated arginine': ('methylation', 'R'),
    'Asymmetric dimethylarginine': ('methylation', 'R'),
    'Omega-N-methylated arginine': ('methylation', 'R'),
    'N6,N6-dimethyllysine': ('methylation', 'K'),
    'N6,N6,N6-trimethyllysine': ('methylation', 'K'),
    'Symmetric dimethylarginine': ('methylation', 'R'),
    'ADP-ribosylarginine': ('ribosylation', 'R'),
    'ADP-ribosylcysteine': ('ribosylation', 'C'),
    'ADP-ribosylasparagine': ('ribosylation', 'N'),
    'PolyADP-ribosyl glutamic acid': ('ribosylation', 'E'),
    'O-acetylserine': ('acetylation', 'S'),
    'O-acetyl-L-serine': ('acetylation', 'S'),
    'N-acetyl-L-alanine': ('acetylation', 'A'),
    'omega-N-methyl-L-arginine': ('methylation', 'R'),
    'symmetric dimethyl-L-arginine': ('methylation', 'R'),
    'N-acetylproline': ('acetylation', 'P'),
    'acetylated': ('acetylation', None),
    'acetylation': ('acetylation', None),
    '(3R)-3-hydroxyaspartate': ('hydroxylation', 'D'),
    '(3R)-3-hydroxyasparagine': ('hydroxylation', 'N'),
    'Tele-methylhistidine': ('methylation', 'H'),
    'N-acetylglutamate': ('acetylation', 'E'),
    'N-acetylaspartate': ('acetylation', 'D'),
    'n6me2lys': ('methylation', 'K'),
    'Phosphohistidine': ('phosphorylation', 'H'),
    'S-farnesyl-L-cysteine': ('farnesylation', 'C'),
    'modified glycine residue': ('modification', 'G'),
    'N-acetyl-L-methionine': ('acetylation', 'M'),
    '4-hydroxy-L-proline': ('hydroxylation', 'P'),
    '4hypro': ('hydroxylation', 'P'),
    '3-hydroxy-L-proline': ('hydroxylation', 'P'),
    '5-hydroxy-L-lysine': ('hydroxylation', 'K'),
    'O-palmitoyl-L-serine': ('palmitoylation', 'S')
    }


def _is_complex(pe):
    """Return True if the physical entity is a complex"""
    return isinstance(pe, bp.Complex)


def _is_protein(pe):
    """Return True if the element is a protein"""
    return isinstance(pe, (bp.Protein, bp.ProteinReference))


def _is_rna(pe):
    """Return True if the element is an RNA"""
    return isinstance(pe, bp.Rna)


def _is_small_molecule(pe):
    """Return True if the element is a small molecule"""
    return isinstance(pe, (bp.SmallMolecule, bp.SmallMoleculeReference))


def _is_physical_entity(pe):
    """Return True if the element is a physical entity"""
    return isinstance(pe, bp.PhysicalEntity)


def _is_simple_physical_entity(pe):
    return isinstance(pe, bp.SimplePhysicalEntity)


def _is_modification(feature):
    return isinstance(feature, bp.ModificationFeature) and \
        get_modification_type(feature) == 'modification'


def _is_activity(feature):
    return isinstance(feature, bp.ModificationFeature) and \
        get_modification_type(feature) in {'activity', 'inactivity'}


def get_modification_type(mf: bp.ModificationFeature):
    if mf.modification_type is None:
        return None

    mf_type_terms = set(mf.modification_type.term)
    if not mf_type_terms:
        return 'modification'

    if mf_type_terms & inactivity_terms:
        return 'inactivity'
    elif mf_type_terms & activity_terms:
        return 'activity'
    else:
        return 'modification'


def _is_reference(bpe):
    """Return True if the element is an entity reference."""
    return isinstance(bpe, bp.EntityReference)


def _is_entity(bpe):
    """Return True if the element is a physical entity."""
    return isinstance(bpe, bp.Entity)


def _is_catalysis(bpe):
    """Return True if the element is Catalysis."""
    return isinstance(bpe, bp.Catalysis)


def _has_members(bpe):
    if _is_reference(bpe):
        members =  bpe.member_entity_reference
    elif _is_entity(bpe):
        members =  bpe.member_physical_entity
    else:
        return False
    if len(members) > 0:
        return True
    else:
        return False


def _listify(lst):
    if not isinstance(lst, collections.Iterable):
        return [lst]
    else:
        return lst


def _list_listify(lst):
    return [l if isinstance(l, collections.Iterable) else [l] for l in lst]


def _get_combinations(lst):
    return itertools.product(*_list_listify(lst))


def _remove_redundant_mods(stmt):
    """Remove redundant Agent states that are modified by statement."""
    # A custom matches function is used so we also match conditions
    # where the polarities are opposite
    def matches(mc1, mc2):
        return mc1.mod_type == mc2.mod_type and \
            mc1.residue == mc2.residue and \
            mc1.position == mc2.position
    stmt_mc = stmt._get_mod_condition()
    stmt.sub = copy.deepcopy(stmt.sub)
    stmt.sub.mods = [mc for mc in stmt.sub.mods if not matches(mc, stmt_mc)]
    return stmt


generic_chebi_ids = {
    '76971',  # E-coli metabolite
    '75771',  # mouse metabolite
    '77746',  # human metabolite
    '27027',  # micronutrient
    '78675',  # fundamental metabolite
    '50860',  # organic molecular entity
}


manual_chebi_map = {
    'H2O': '15377',
    'phosphate': '18367',
    'H+': '15378',
    'O2': '15379'
}


@lru_cache(maxsize=5000)
def get_specific_chebi_id(chebi_ids, name):
    # NOTE: this function is mainly factored out to be able to use cacheing, it
    # requires a frozenset as input to work.

    # First, if we have a manual override, we just do that
    manual_id = manual_chebi_map.get(name)
    if manual_id:
        return manual_id

    # The first thing we do is eliminate the secondary IDs by mapping them to
    # primaries
    primary_ids = {chebi_client.get_primary_id(cid)
                   for cid in chebi_ids}
    # Occasinally, invalid ChEBI IDs are given that don't have corresponding
    # primary IDs, which we can filter out
    primary_ids = {pi for pi in primary_ids if pi is not None}
    # We then get rid of generic IDs which are never useful for grounding
    non_generic_ids = primary_ids - generic_chebi_ids

    # We then try name-based grounding to see if any of the names in the list
    # match the name of the entity well enough
    grounding_names = [chebi_client.get_chebi_name_from_id(p) for p in
                       non_generic_ids]
    for grounding_name, grounding_id in zip(grounding_names, non_generic_ids):
        if grounding_name and (name.lower() == grounding_name.lower()):
            return grounding_id

    # If we still have no best grounding, we try to distill the IDs down to
    # the most specific one based on the hierarchy
    specific_chebi_id = chebi_client.get_specific_id(non_generic_ids)
    return specific_chebi_id


activity_terms = {
    'active',
    'residue modification, active',
}


inactivity_terms = {
    'inactive',
    'residue modification, inactive',
}


