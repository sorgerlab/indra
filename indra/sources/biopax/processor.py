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
from indra.ontology.standardize import standardize_name_db_refs

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
    def __init__(self, model, use_conversion_level_evidence=True):
        self.model = model
        self.statements = []
        self._mod_conditions = {}
        self._activity_conditions = {}
        self._agents = {}
        self.use_conversion_level_evidence = use_conversion_level_evidence

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

    def feature_delta(self, from_pe: bp.PhysicalEntity,
                      to_pe: bp.PhysicalEntity):
        """Return gained and lost modifications and any activity change."""
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

    def _extract_features(self):
        """Pre-extract features before processing statements."""
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

    @staticmethod
    def find_matching_left_right(conversion: bp.Conversion):
        """Find matching entities on the left and right of a conversion."""
        left_simple = []
        for pe in conversion.left:
            left_simple += expand_complex(pe)
        right_simple = []
        for pe in conversion.right:
            right_simple += expand_complex(pe)

        return BiopaxProcessor.find_matching_entities(left_simple, right_simple)

    @staticmethod
    def find_matching_entities(left_simple, right_simple):
        """Find matching entities between two lists of simple entities."""
        matches = []
        for inp, outp in itertools.product(left_simple, right_simple):
            inp_type = infer_pe_type(inp)
            outp_type = infer_pe_type(outp)
            if inp_type != outp_type:
                continue
            elif inp_type == 'family':
                input_ers = {mpe.entity_reference.uid
                             for mpe in expand_family(inp)
                             if isinstance(mpe, bp.SimplePhysicalEntity)}
                output_ers = {mpe.entity_reference.uid
                              for mpe in expand_family(outp)
                              if isinstance(mpe, bp.SimplePhysicalEntity)}
                if input_ers == output_ers:
                    matches.append((inp, outp))
            # Sometimes we get a "raw" PhysicalEntity here that doesn't
            # actually have an entity reference which is required here so
            # we skip those.
            elif not isinstance(inp, (bp.SimplePhysicalEntity, bp.Complex)) or \
                    not isinstance(outp, (bp.SimplePhysicalEntity, bp.Complex)):
                continue
            elif inp_type == 'complex_named':
                if inp.uid == outp.uid:
                    matches.append((inp, outp))
            elif inp_type == 'complex_family':
                inp_members = flatten([expand_complex(m)
                                      for m in expand_family(inp)])
                outp_members = flatten([expand_complex(m)
                                        for m in expand_family(outp)])
                matches += BiopaxProcessor.find_matching_entities(inp_members,
                                                                  outp_members)
            elif inp.entity_reference.uid == outp.entity_reference.uid:
                matches.append((inp, outp))
        return matches

    def _control_conversion_iter(self, conversion_type, controller_logic):
        """An iterator over controlled conversions in the model."""
        for control in self.model.get_objects_by_type(bp.Control):
            conversion = control.controlled
            # Sometimes there is nothing being controlled, we skip
            # those cases. We also want to skip things like Modulation
            # that don't convert anything.
            if not isinstance(conversion, conversion_type):
                continue
            ev = self._get_evidence(control)
            for controller_pe in control.controller:
                # We skip e.g., Pathway controllers
                if not isinstance(controller_pe, bp.PhysicalEntity):
                    continue
                if controller_logic == 'primary':
                    primary_controller = \
                        self._get_primary_controller(controller_pe)
                else:
                    primary_controller = \
                        self._get_all_protein_controllers(controller_pe)
                if not primary_controller:
                    continue
                primary_controller_agents = []
                for pc in _listify(primary_controller):
                    primary_controller_agents += \
                        _listify(self._get_agents_from_entity(pc))
                for primary_controller_agent in primary_controller_agents:
                    yield primary_controller_agent, ev, control, conversion

    def _conversion_no_control_iter(self):
        """An iterator over conversions irrespective of control in the model."""
        for conversion in self.model.get_objects_by_type(bp.Conversion):
            ev = self._get_evidence(conversion)
            for inp, outp in self.find_matching_left_right(conversion):
                for inp_simple, outp_simple in \
                        zip(expand_family(inp),
                            expand_family(outp)):
                    gained_mods, lost_mods, activity_change = \
                        self.feature_delta(inp_simple, outp_simple)
                    inp_agents = \
                        self._get_agents_from_singular_entity(inp_simple)
                    for inp_agent in inp_agents:
                        yield inp_agent, gained_mods, lost_mods, \
                            activity_change, ev

    def _conversion_state_iter(self):
        """An iterator over state changed in controlled conversions
        in the model."""
        for primary_controller_agent, ev, control, conversion in \
                self._control_conversion_iter(bp.Conversion, 'primary'):
            for inp, outp in self.find_matching_left_right(conversion):
                # There is sometimes activity change at the family level
                # which we need to capture
                _, _, overall_activity_change = self.feature_delta(inp, outp)
                for inp_simple, outp_simple in \
                        zip(expand_family(inp),
                            expand_family(outp)):
                    gained_mods, lost_mods, activity_change = \
                        self.feature_delta(inp_simple, outp_simple)
                    activity_change = activity_change if activity_change else \
                        overall_activity_change
                    inp_agents = \
                        self._get_agents_from_singular_entity(inp_simple)
                    for inp_agent in inp_agents:
                        yield primary_controller_agent, inp_agent, \
                              gained_mods, lost_mods, activity_change, ev

    def get_modifications(self):
        """Extract INDRA Modification Statements from the BioPAX model."""
        for enz, sub, gained_mods, lost_mods, \
                activity_change, ev in self._conversion_state_iter():
            for mods, is_gain in ((gained_mods, True), (lost_mods, False)):
                for mod in mods:
                    stmt_class = modtype_to_modclass[mod.mod_type]
                    if not is_gain:
                        stmt_class = modclass_to_inverse[stmt_class]
                    stmt = stmt_class(enz, sub, mod.residue,
                                      mod.position, evidence=ev)
                    stmt = _remove_redundant_mods(stmt)
                    self.statements.append(stmt)

    def get_regulate_activities(self):
        """Get Activation/Inhibition INDRA Statements from the BioPAX model."""
        for subj, obj, gained_mods, lost_mods, \
                activity_change, ev in self._conversion_state_iter():
            # We don't want to have gained or lost modification features
            if not activity_change or (gained_mods or lost_mods):
                continue
            stmt_class = Activation if activity_change == 'active' \
                else Inhibition
            stmt = stmt_class(subj, obj, 'activity',
                              evidence=ev)
            self.statements.append(stmt)

    def get_activity_modification(self):
        """Extract INDRA ActiveForm statements from the BioPAX model."""
        for agent, gained_mods, lost_mods, activity_change, ev in \
                self._conversion_no_control_iter():
            # We have to have both a modification change and an activity
            # change
            if not (gained_mods or lost_mods) or not activity_change:
                continue
            is_active = (activity_change == 'active')
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
        for subj, ev, control, conversion in \
                self._control_conversion_iter(bp.TemplateReaction, 'all'):
            stmt_type = IncreaseAmount if control.control_type == 'ACTIVATION' \
                else DecreaseAmount
            for product in conversion.product:
                product_agents = self._get_agents_from_entity(product)
                for obj in _listify(product_agents):
                    stmt = stmt_type(subj, obj, evidence=ev)
                    self.statements.append(stmt)

    def get_conversions(self):
        """Extract Conversion INDRA Statements from the BioPAX model."""
        for subj, ev, control, conversion in \
                self._control_conversion_iter(bp.Conversion, 'primary'):
            # We only extract conversions for small molecules
            if not all(_is_small_molecule(pe)
                       for pe in (conversion.left + conversion.right)):
                continue
            # Since we don't extract location, this produces conversions where
            # the input and output is the same
            if isinstance(conversion, bp.Transport):
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
            st = Conversion(subj, obj_from, obj_to, evidence=ev)
            self.statements.append(st)

    @staticmethod
    def find_gdp_gtp_complex(cplxes):
        for cplx in cplxes:
            members = expand_complex(cplx)
            if not members or len(members) != 2:
                continue
            gdp_gtp_idx = None
            ras_member = None
            for idx, member in enumerate(members):
                if _is_small_molecule(member) and \
                        member.display_name in {'GDP', 'GTP'}:
                    gdp_gtp_idx = idx
                elif _is_protein(member):
                    ras_member = member
            if gdp_gtp_idx is None or ras_member is None:
                continue
            return ras_member, members[gdp_gtp_idx].display_name
        return None, None

    def get_gap_gef(self):
        """Extract Gap and Gef INDRA Statements."""
        for gap_gef, ev, control, conversion in \
                self._control_conversion_iter(bp.Conversion, 'primary'):
            assert isinstance(gap_gef, Agent)
            # We have to make sure that we don't pick up chemicals here
            if not set(gap_gef.db_refs.keys()) & {'HGNC', 'UP'}:
                continue
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

            ras_agents = self._get_agents_from_entity(left_ras)
            for ras in _listify(ras_agents):
                st = stmt_type(gap_gef, ras, evidence=ev)
                self.statements.append(st)

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
        # If it's not a real complex, just return the physical entity
        # as is
        if infer_pe_type(controller_pe) != 'complex':
            return controller_pe

        # Identifying the "real" enzyme in a complex may not always be
        # possible.
        # One heuristic here could be to find the member which is
        # active and if it is the only active member then
        # set this as the enzyme to which all other members of the
        # complex are bound.
        # Separate out protein and non-protein members
        protein_members = [p for p in expand_complex(controller_pe)
                           if _is_protein(p)]
        # If there is only one protein member, we can assume that
        # it is the enzyme, and everything else is just bound
        # to it.
        if len(protein_members) == 1:
            return protein_members[0]
        return None

    def _get_all_protein_controllers(self, controller_pe):
        # If it's not a real complex, just return the physical entity
        # as is
        if infer_pe_type(controller_pe) != 'complex':
            return controller_pe

        protein_members = [p for p in controller_pe.component
                           if _is_protein(p)]
        return protein_members

    def _get_agents_from_singular_entity(self, bpe: bp.PhysicalEntity):
        """This is for extracting one or more Agents from a PhysicalEntity
        which doesn't have member_physical_entities."""
        try:
            return copy.deepcopy(self._agents[bpe.uid])
        except KeyError:
            pass

        mcs = BiopaxProcessor._get_entity_mods(bpe) if _is_protein(bpe) else []
        name = bpe.display_name
        agents = []

        # We first get processed xrefs
        xrefs = BiopaxProcessor._get_processed_xrefs(bpe)

        # We now need to harmonize UP and HGNC
        # Case 1. Multiple genes coding for one protein
        nhgnc_ids = len(xrefs.get('HGNC', {}))
        nup_ids = len(xrefs.get('UP', {}))
        # One protein coded by many genes
        if nhgnc_ids > 1 and nup_ids == 1:
            for hgnc_id in xrefs['HGNC']:
                standard_name, db_refs = \
                    standardize_name_db_refs({'HGNC': hgnc_id})
                if standard_name:
                    name = standard_name
                agents.append(Agent(name, db_refs=db_refs, mods=mcs))
        # One gene coding for many proteins
        elif nhgnc_ids == 1 and nup_ids > 1:
            for up_id in xrefs['UP']:
                standard_name, db_refs = \
                    standardize_name_db_refs({'UP': up_id})
                if standard_name:
                    name = standard_name
                agents.append(Agent(name, db_refs=db_refs, mods=mcs))
        # This is secretly a family, i.e., we have more than one
        # gene/protein IDs and so we can go by one of the ID sets and
        # standardize from there
        elif nhgnc_ids > 1 and nhgnc_ids == nup_ids:
            for up_id in xrefs['UP']:
                standard_name, db_refs = \
                    standardize_name_db_refs({'UP': up_id})
                if standard_name:
                    name = standard_name
                agents.append(Agent(name, db_refs=db_refs, mods=mcs))
        # Otherwise it's just a regular Agent
        else:
            standard_name, db_refs = \
                standardize_name_db_refs(clean_up_xrefs(xrefs))
            if standard_name:
                name = standard_name
            agents.append(Agent(name, db_refs=db_refs, mods=mcs))
        return agents

    def _get_agents_from_entity(self, bpe: bp.PhysicalEntity):
        # If the entity has members (like a protein family),
        # we iterate over them
        if bpe.member_physical_entity:
            agents = []
            for m in bpe.member_physical_entity:
                member_agents = self._get_agents_from_entity(m)
                agents += member_agents
            return agents

        # If it is a single entity, we get its name and database
        # references
        return self._get_agents_from_singular_entity(bpe)

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

    def _get_evidence(self, bpe: bp.PhysicalEntity):
        citations = BiopaxProcessor._get_citations(bpe)
        if self.use_conversion_level_evidence and hasattr(bpe, 'controller'):
            citations += BiopaxProcessor._get_citations(bpe.controlled)
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
        source_id = 'http://pathwaycommons.org/pc12/%s' % bpe.uid if \
            not bpe.uid.startswith('http') else bpe.uid
        ev = [Evidence(source_api='biopax', pmid=cit,
                       source_id=source_id, epistemics=epi,
                       annotations=annotations)
              for cit in citations]
        return ev

    @staticmethod
    def _get_citations(bpe: bp.PhysicalEntity):
        refs = []
        xrefs = bpe.xref + flatten([ev.xref for ev in bpe.evidence])
        for xr in xrefs:
            db_name = xr.db
            if db_name is not None and db_name.upper() == 'PUBMED':
                refs.append(xr.id)
        # TODO: handle non-pubmed evidence
        return refs

    @staticmethod
    def _get_processed_xrefs(bpe: bp.PhysicalEntity):
        entref = BiopaxProcessor._get_entref(bpe)
        if not entref:
            return {}

        primary_ns, primary_id = \
            BiopaxProcessor._get_reference_primary_id(entref)

        from collections import defaultdict
        xrefs = defaultdict(set)
        if primary_ns and primary_id:
            xrefs[primary_ns].add(primary_id)

        for xref in entref.xref:
            if not xref.db:
                continue
            xref_db_ns = xref_ns_map.get(xref.db.lower())
            if not xref_db_ns:
                continue
            xrefs[xref_db_ns].add(xref.id)

        xrefs = dict(xrefs)

        # We now sanitize certain key name spaces
        if 'UP' in xrefs:
            up_ids = sanitize_up_ids(xrefs['UP'])
            if up_ids:
                xrefs['UP'] = up_ids
            else:
                xrefs.pop('UP')
        if 'HGNC' in xrefs or 'HGNC.SYMBOL' in xrefs:
            hgnc_ids = xrefs.get('HGNC', set()) | \
                xrefs.get('HGNC.SYMBOL', set())
            hgnc_ids = sanitize_hgnc_ids(hgnc_ids)
            if hgnc_ids:
                xrefs['HGNC'] = hgnc_ids
            else:
                xrefs.pop('HGNC', None)
            xrefs.pop('HGNC.SYMBOL', None)
        if 'CHEBI' in xrefs:
            chebi_id = sanitize_chebi_ids(xrefs['CHEBI'], bpe.display_name)
            if chebi_id:
                xrefs['CHEBI'] = chebi_id
            else:
                xrefs.pop('CHEBI')

        return xrefs

    @staticmethod
    def _get_reference_primary_id(entref: bp.EntityReference):
        # In practice, it appears that only UniProt and ChEBI appear in this
        # form.
        match = re.match('http://identifiers.org/([^/]+)/(.+)$',
                         entref.uid)
        if match:
            ident_ns, ident_id = match.groups()
            if ident_ns == 'uniprot':
                primary_ns, primary_id = 'UP', ident_id
            elif ident_ns == 'chebi':
                primary_ns, primary_id = 'CHEBI', ident_id
            elif ident_ns == 'pubchem.compound':
                primary_ns, primary_id = 'PUBCHEM', ident_id
            elif ident_ns == 'pubchem.substance':
                primary_ns, primary_id = 'PUBCHEM.SUBSTANCE', ident_id
            else:
                logger.warning('Unhandled identifiers namespace: %s' %
                               ident_ns)
                primary_ns, primary_id = None, None
        else:
            primary_ns, primary_id = None, None

        return primary_ns, primary_id

    @staticmethod
    def _get_entref(bpe: bp.PhysicalEntity):
        """Returns the entity reference of an entity if it exists or
        return the entity reference that was passed in as argument."""
        if isinstance(bpe, bp.SimplePhysicalEntity):
            return bpe.entity_reference
        return None

    def get_coverage(self):
        uids = set()
        for uid, obj in self.model.objects.items():
            if isinstance(obj, (bp.Catalysis, bp.TemplateReactionRegulation)):
                uids.add(uid)
        stmt_uids = set()
        for stmt in self.statements:
            for ev in stmt.evidence:
                stmt_uids.add(ev.source_id)

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


def _is_small_molecule(pe):
    """Return True if the element is a small molecule"""
    return isinstance(pe, (bp.SmallMolecule, bp.SmallMoleculeReference))


def _is_modification(feature):
    return isinstance(feature, bp.ModificationFeature) and \
        get_modification_type(feature) == 'modification'


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


def _is_entity(bpe):
    """Return True if the element is a physical entity."""
    return isinstance(bpe, bp.Entity)


def _listify(lst):
    if not isinstance(lst, collections.Iterable):
        return [lst]
    else:
        return lst


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


def expand_complex(pe: bp.PhysicalEntity):
    simple_pes = []
    if _is_complex(pe) and pe.component:
        for component_pe in pe.component:
            simple_pes += expand_complex(component_pe)
    else:
        simple_pes.append(pe)
    return simple_pes


def expand_family(pe: bp.PhysicalEntity):
    if pe.member_physical_entity:
        return pe.member_physical_entity
    else:
        return [pe]


def infer_pe_type(pe: bp.PhysicalEntity):
    if isinstance(pe, bp.Complex):
        if pe.component:
            return 'complex'
        elif pe.member_physical_entity:
            return 'complex_family'
        else:
            return 'complex_named'
    else:
        if pe.member_physical_entity:
            return 'family'
        else:
            return 'single'


generic_chebi_ids = {
    'CHEBI:76971',  # E-coli metabolite
    'CHEBI:75771',  # mouse metabolite
    'CHEBI:77746',  # human metabolite
    'CHEBI:27027',  # micronutrient
    'CHEBI:78675',  # fundamental metabolite
    'CHEBI:50860',  # organic molecular entity
}


manual_chebi_map = {
    'H2O': 'CHEBI:15377',
    'phosphate': 'CHEBI:18367',
    'H+': 'CHEBI:15378',
    'O2': 'CHEBI:15379'
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


def sanitize_up_ids(up_ids):
    # First, we map any secondary IDs to primary IDs
    up_ids = {uniprot_client.get_primary_id(up_id)
              for up_id in up_ids}
    # We filter out IDs that are actually mnemonics, these are just mixed
    # in without any differentiation from other IDs
    up_ids = {up_id for up_id in up_ids if '_' not in up_id}
    # TODO: should we do anything about isoforms?
    # We separate out specific sets of IDs
    human_ids = [up_id for up_id in up_ids
                 if uniprot_client.is_human(up_id)]
    reviewed_non_human_ids = [
        up_id for up_id in up_ids
        if not uniprot_client.is_human(up_id)
        # get_mnemonic is just a quick way to see if we have this entry
        and uniprot_client.get_mnemonic(up_id, web_fallback=False)]
    if human_ids:
        return human_ids
    elif reviewed_non_human_ids:
        return reviewed_non_human_ids
    else:
        return []


def sanitize_chebi_ids(chebi_ids, name):
    chebi_ids = {chebi_id if chebi_id.startswith('CHEBI:')
                 else 'CHEBI:%s' % chebi_id for chebi_id in chebi_ids}
    chebi_ids = {chebi_client.get_primary_id(chebi_id)
                 for chebi_id in chebi_ids}
    if len(chebi_ids) == 1:
        return list(chebi_ids)
    specific_chebi_id = get_specific_chebi_id(frozenset(chebi_ids),
                                              name)
    return specific_chebi_id


def sanitize_hgnc_ids(raw_hgnc_ids):
    # First we get a list of primary IDs
    hgnc_ids = set()
    for raw_hgnc_id in raw_hgnc_ids:
        # Check if it's an ID first
        m1 = re.match('([0-9]+)', raw_hgnc_id)
        m2 = re.match('hgnc:([0-9]+)', raw_hgnc_id.lower())
        if m1:
            hgnc_id = str(m1.groups()[0])
            hgnc_ids.add(hgnc_id)
        elif m2:
            hgnc_id = str(m2.groups()[0])
            hgnc_ids.add(hgnc_id)
        # If not, we assume it's a symbol
        else:
            hgnc_id = hgnc_client.get_current_hgnc_id(raw_hgnc_id)
            if isinstance(hgnc_id, list):
                hgnc_ids |= set(hgnc_id)
            elif hgnc_id:
                hgnc_ids.add(hgnc_id)

    return list(hgnc_ids)


def clean_up_xrefs(xrefs):
    db_refs = {}
    for k, v in xrefs.items():
        if isinstance(v, (list, set)):
            if len(v) == 1:
                db_refs[k] = list(v)[0]
        else:
            db_refs[k] = v
    return db_refs


activity_terms = {
    'active',
    'residue modification, active',
}


inactivity_terms = {
    'inactive',
    'residue modification, inactive',
}


xref_ns_map = {
    'chebi': 'CHEBI',
    'uniprot': 'UP',
    'uniprot isoform': 'UP',
    'uniprot knowledgebase': 'UP',
    'uniprotkb': 'UP',
    'ncbi gene': 'EGID',
    'hgnc symbol': 'HGNC.SYMBOL',
    'hgnc': 'HGNC',
    'mesh': 'MESH',
    'drugbank': 'DRUGBANK',
    'pubchem-compound': 'PUBCHEM',
    'chembl compound': 'CHEMBL',
    'pubchem': 'PUBCHEM',
    'cas': 'CAS',
    'hmdb': 'HMDB',
    'gene ontology': 'GO',
    'interpro': 'IP',
    'mirbase': 'MIRBASE',
    'mirbase mature sequence': 'MIRBASEM',
    'hugo gene nomenclature committee (hgnc)': 'HGNC',
    'ensembl': 'ENSEMBL',
    'taxonomy': 'TAXONOMY',
}
