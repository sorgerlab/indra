from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import sys
import logging
import itertools
import collections
try:
    from functools import lru_cache
except ImportError:
    from functools32 import lru_cache

from indra.java_vm import autoclass, JavaException, cast
from indra.databases import hgnc_client, uniprot_client, chebi_client
from indra.statements import *
from . import pathway_commons_client as pcc
from indra.util import decode_obj

logger = logging.getLogger('biopax')

# TODO:
# - Extract cellularLocation from each PhysicalEntity
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

    def print_statements(self):
        """Print all INDRA Statements collected by the processors."""
        for i, stmt in enumerate(self.statements):
            print("%s: %s" % (i, stmt))

    def save_model(self, file_name=None):
        """Save the BioPAX model object in an OWL file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the OWL file to save the model in.
        """
        if file_name is None:
            logger.error('Missing file name')
            return
        pcc.model_to_owl(self.model, file_name)

    def get_complexes(self):
        """Extract INDRA Complex Statements from the BioPAX model.

        This method searches for org.biopax.paxtools.model.level3.Complex
        objects which represent molecular complexes. It doesn't reuse
        BioPAX Pattern's org.biopax.paxtools.pattern.PatternBox.inComplexWith
        query since that retrieves pairs of complex members rather than
        the full complex.
        """
        for obj in self.model.getObjects().toArray():
            bpe = _cast_biopax_element(obj)
            if not _is_complex(bpe):
                continue
            ev = self._get_evidence(bpe)

            members = self._get_complex_members(bpe)
            if members is not None:
                if len(members) > 10:
                    logger.debug('Skipping complex with more than 10 members.')
                    continue
                complexes = _get_combinations(members)
                for c in complexes:
                    self.statements.append(decode_obj(Complex(c, ev),
                                                      encoding='utf-8'))

    def get_modifications(self):
        """Extract INDRA Modification Statements from the BioPAX model.

        To extract Modifications, this method reuses the structure of
        BioPAX Pattern's
        org.biopax.paxtools.pattern.PatternBox.constrolsStateChange pattern
        with additional constraints to specify the type of state change
        occurring (phosphorylation, deubiquitination, etc.).
        """
        for modtype, modclass in modtype_to_modclass.items():
            # TODO: we could possibly try to also extract generic
            # modifications here
            if modtype == 'modification':
                continue
            stmts = self._get_generic_modification(modclass)
            self.statements += stmts

    def get_activity_modification(self):
        """Extract INDRA ActiveForm statements from the BioPAX model.

        This method extracts ActiveForm Statements that are due to
        protein modifications. This method reuses the structure of
        BioPAX Pattern's
        org.biopax.paxtools.pattern.PatternBox.constrolsStateChange pattern
        with additional constraints to specify the gain or loss of a
        modification occurring (phosphorylation, deubiquitination, etc.)
        and the gain or loss of activity due to the modification state
        change.
        """
        mod_filter = 'residue modification, active'
        for is_active in [True, False]:
            p = self._construct_modification_pattern()
            rel = mcct.GAIN if is_active else mcct.LOSS
            p.add(mcc(rel, mod_filter),
                  "input simple PE", "output simple PE")

            s = _bpp('Searcher')
            res = s.searchPlain(self.model, p)
            res_array = [_match_to_array(m) for m in res.toArray()]

            for r in res_array:
                reaction = r[p.indexOf('Conversion')]
                activity = 'activity'
                input_spe = r[p.indexOf('input simple PE')]
                output_spe = r[p.indexOf('output simple PE')]

                # Get the modifications
                mod_in = \
                    BiopaxProcessor._get_entity_mods(input_spe)
                mod_out = \
                    BiopaxProcessor._get_entity_mods(output_spe)

                mod_shared = _get_mod_intersection(mod_in, mod_out)
                gained_mods = _get_mod_difference(mod_out, mod_in)

                # Here we get the evidence for the BiochemicalReaction
                ev = self._get_evidence(reaction)

                agents = self._get_agents_from_entity(output_spe)
                for agent in _listify(agents):
                    static_mods = _get_mod_difference(agent.mods,
                                                      gained_mods)
                    # NOTE: with the ActiveForm representation we cannot
                    # separate static_mods and gained_mods. We assume here
                    # that the static_mods are inconsequential and therefore
                    # are not mentioned as an Agent condition, following
                    # don't care don't write semantics. Therefore only the
                    # gained_mods are listed in the ActiveForm as Agent
                    # conditions.
                    if gained_mods:
                        agent.mods = gained_mods
                        stmt = ActiveForm(agent, activity, is_active,
                                          evidence=ev)
                        self.statements.append(decode_obj(stmt,
                                                          encoding='utf-8'))

    def get_regulate_activities(self):
        """Get Activation/Inhibition INDRA Statements from the BioPAX model.

        This method extracts Activation/Inhibition Statements and reuses the
        structure of BioPAX Pattern's
        org.biopax.paxtools.pattern.PatternBox.constrolsStateChange pattern
        with additional constraints to specify the gain or loss of
        activity state but assuring that the activity change is not due to
        a modification state change (which are extracted by get_modifications
        and get_activity_modification).
        """
        mcc = _bpp('constraint.ModificationChangeConstraint')
        mcct = _bpp('constraint.ModificationChangeConstraint$Type')
        mod_filter = 'residue modification, active'
        # Start with a generic modification pattern
        p = BiopaxProcessor._construct_modification_pattern()
        stmts = []
        for act_class, gain_loss in zip([Activation, Inhibition],
                                        [mcct.GAIN, mcct.LOSS]):
            p.add(mcc(gain_loss, mod_filter),
                      "input simple PE", "output simple PE")
            s = _bpp('Searcher')
            res = s.searchPlain(self.model, p)
            res_array = [_match_to_array(m) for m in res.toArray()]
            for r in res_array:
                controller_pe = r[p.indexOf('controller PE')]
                input_pe = r[p.indexOf('input PE')]
                input_spe = r[p.indexOf('input simple PE')]
                output_spe = r[p.indexOf('output simple PE')]
                reaction = r[p.indexOf('Conversion')]
                control = r[p.indexOf('Control')]

                if not _is_catalysis(control):
                    continue
                cat_dir = control.getCatalysisDirection()
                if cat_dir is not None and cat_dir.name() != 'LEFT_TO_RIGHT':
                    logger.debug('Unexpected catalysis direction: %s.' % \
                        control.getCatalysisDirection())
                    continue

                subjs = BiopaxProcessor._get_primary_controller(controller_pe)
                if not subjs:
                    continue

                '''
                if _is_complex(input_pe):
                    # TODO: It is possible to find which member of the complex
                    # is actually activated. That member will be the substrate
                    # and all other members of the complex will be bound to it.
                    logger.info('Cannot handle complex subjects.')
                    continue
                '''
                objs = BiopaxProcessor._get_agents_from_entity(input_spe,
                                                               expand_pe=False)

                ev = self._get_evidence(control)
                for subj, obj in itertools.product(_listify(subjs),
                                                   _listify(objs)):
                    # Get the modifications
                    mod_in = \
                        BiopaxProcessor._get_entity_mods(input_spe)
                    mod_out = \
                        BiopaxProcessor._get_entity_mods(output_spe)

                    # We assume if modifications change then this is not really
                    # a pure activation event
                    gained_mods = _get_mod_difference(mod_out, mod_in)
                    lost_mods = _get_mod_difference(mod_in, mod_out)
                    if gained_mods or lost_mods:
                        continue

                    stmt = act_class(subj, obj, 'activity', evidence=ev)
                    self.statements.append(decode_obj(stmt, encoding='utf-8'))


    def get_regulate_amounts(self):
        """Extract INDRA RegulateAmount Statements from the BioPAX model.

        This method extracts IncreaseAmount/DecreaseAmount Statements from
        the BioPAX model. It fully reuses BioPAX Pattern's
        org.biopax.paxtools.pattern.PatternBox.controlsExpressionWithTemplateReac
        pattern to find TemplateReactions which control the expression of
        a protein.
        """
        p = pb.controlsExpressionWithTemplateReac()
        s = _bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [_match_to_array(m) for m in res.toArray()]
        stmts = []
        for res in res_array:
            # FIXME: for some reason labels are not accessible
            # for these queries. It would be more reliable
            # to get results by label instead of index.
            '''
            controller_er = res[p.indexOf('controller ER')]
            generic_controller_er = res[p.indexOf('generic controller ER')]
            controller_simple_pe = res[p.indexOf('controller simple PE')]
            controller_pe = res[p.indexOf('controller PE')]
            control = res[p.indexOf('Control')]
            conversion = res[p.indexOf('Conversion')]
            input_pe = res[p.indexOf('input PE')]
            input_simple_pe = res[p.indexOf('input simple PE')]
            changed_generic_er = res[p.indexOf('changed generic ER')]
            output_pe = res[p.indexOf('output PE')]
            output_simple_pe = res[p.indexOf('output simple PE')]
            changed_er = res[p.indexOf('changed ER')]
            '''
            # TODO: here, res[3] is the complex physical entity
            # for instance http://pathwaycommons.org/pc2/
            # Complex_43c6b8330562c1b411d21e9d1185bae9
            # consists of 3 components: JUN, FOS and NFAT
            # where NFAT further contains 3 member physical entities.
            #
            # However, res[2] iterates over all 5 member physical entities
            # of the complex which doesn't represent the underlying
            # structure faithfully. It would be better to use res[3]
            # (the complex itself) and look at components and then
            # members. However, then, it would not be clear how to
            # construct an INDRA Agent for the controller.
            controller = self._get_agents_from_entity(res[2])
            controlled_pe = res[6]
            controlled = self._get_agents_from_entity(controlled_pe)

            conversion = res[5]
            direction = conversion.getTemplateDirection()
            if direction is not None:
                direction = direction.name()
                if direction != 'FORWARD':
                    logger.warning('Unhandled conversion direction %s' %
                                   direction)
                    continue
            # Sometimes interaction type is annotated as
            # term=='TRANSCRIPTION'. Other times this is not
            # annotated.
            int_type = conversion.getInteractionType().toArray()
            if int_type:
                for it in int_type:
                    for term in it.getTerm().toArray():
                        pass
            control = res[4]
            control_type = control.getControlType()
            if control_type:
                control_type = control_type.name()
            ev = self._get_evidence(control)
            for subj, obj in itertools.product(_listify(controller),
                                               _listify(controlled)):
                subj_act = ActivityCondition('transcription', True)
                subj.activity = subj_act
                if control_type == 'ACTIVATION':
                    st = IncreaseAmount(subj, obj, evidence=ev)
                elif control_type == 'INHIBITION':
                    st = DecreaseAmount(subj, obj, evidence=ev)
                else:
                    logger.warning('Unhandled control type %s' % control_type)
                    continue
                st_dec = decode_obj(st, encoding='utf-8')
                self.statements.append(st_dec)

    def get_conversions(self):
        """Extract Conversion INDRA Statements from the BioPAX model.

        This method uses a custom BioPAX Pattern
        (one that is not implemented PatternBox) to query for
        BiochemicalReactions whose left and right hand sides are collections
        of SmallMolecules. This pattern thereby extracts metabolic
        conversions as well as signaling processes via small molecules
        (e.g. lipid phosphorylation or cleavage).
        """
        # NOTE: This pattern gets all reactions in which a protein is the
        # controller and chemicals are converted. But with this pattern only
        # a single chemical is extracted from each side. This can be misleading
        # since we want to capture all inputs and all outputs of the
        # conversion. So we need to step back to the conversion itself and
        # enumerate all inputs/outputs, make sure they constitute the kind
        # of conversion we can capture here and then extract as a Conversion
        # Statement. Another issue here is that the same reaction will be
        # extracted multiple times if there is more then one input or output.
        # Therefore we need to cache the ID of the reactions that have already
        # been handled.
        p = _bpp('Pattern')(_bpimpl('PhysicalEntity')().getModelInterface(),
                            'controller PE')
        # Getting the control itself
        p.add(cb.peToControl(), "controller PE", "Control")
        # Make sure the controller is a protein
        # TODO: possibly allow Complex too
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "controller PE")
        # Link the control to the conversion that it controls
        p.add(cb.controlToConv(), "Control", "Conversion")
        # Make sure this is a BiochemicalRection (as opposed to, for instance,
        # ComplexAssembly)
        p.add(tp(_bpimpl('BiochemicalReaction')().getModelInterface()),
                         "Conversion")
        # The controller shouldn't be a participant of the conversion
        p.add(_bpp('constraint.NOT')(cb.participant()),
              "Conversion", "controller PE")
        # Get the input participant of the conversion
        p.add(pt(rt.INPUT, True), "Control", "Conversion", "input PE")
        # Link to the other side of the conversion
        p.add(cs(cst.OTHER_SIDE), "input PE", "Conversion", "output PE")
        # Make sure the two sides are not the same
        p.add(_bpp('constraint.Equality')(False), "input PE", "output PE")
        # Make sure the input/output is a chemical
        p.add(tp(_bpimpl('SmallMolecule')().getModelInterface()), "input PE")
        p.add(tp(_bpimpl('SmallMolecule')().getModelInterface()), "output PE")

        s = _bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [_match_to_array(m) for m in res.toArray()]
        stmts = []
        reaction_extracted = set()
        for r in res_array:
            controller_pe = r[p.indexOf('controller PE')]
            reaction = r[p.indexOf('Conversion')]
            control = r[p.indexOf('Control')]
            input_pe = r[p.indexOf('input PE')]
            output_pe = r[p.indexOf('output PE')]
            if control.getUri() in reaction_extracted:
                continue
            # Get controller
            subj_list = self._get_agents_from_entity(controller_pe)
            # Get inputs and outputs
            left = reaction.getLeft().toArray()
            right = reaction.getRight().toArray()
            # Skip this if not all participants are chemicals
            if any([not _is_small_molecule(e) for e in left]):
                continue
            if any([not _is_small_molecule(e) for e in right]):
                continue
            obj_left = []
            obj_right = []
            for participant in left:
                agent = self._get_agents_from_entity(participant)
                obj_left.append(agent)
            for participant in right:
                agent = self._get_agents_from_entity(participant)
                obj_right.append(agent)
            ev = self._get_evidence(control)
            for subj in _listify(subj_list):
                st = Conversion(subj, obj_left, obj_right, evidence=ev)
                st_dec = decode_obj(st, encoding='utf-8')
                self.statements.append(st_dec)
            reaction_extracted.add(control.getUri())

    def _gef_gap_base(self):
        # The following constraints were pieced together based on the
        # following two higher level constrains: pb.controlsStateChange(),
        # pb.controlsPhosphorylation().
        p = _bpp('Pattern')(_bpimpl('PhysicalEntity')().getModelInterface(),
                            'controller PE')
        # Getting the control itself
        p.add(cb.peToControl(), "controller PE", "Control")
        # Link the control to the conversion that it controls
        p.add(cb.controlToConv(), "Control", "Conversion")
        # The controller shouldn't be a participant of the conversion
        p.add(_bpp('constraint.NOT')(cb.participant()),
              "Conversion", "controller PE")
        # Get the input participant of the conversion
        p.add(pt(rt.INPUT, True), "Control", "Conversion", "input PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "input PE", "input simple PE")
        # Link to ER
        p.add(cb.peToER(), "input simple PE", "input simple ER")
        # Make sure the participant is a protein
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "input simple PE")
        # Make sure the controller is a protein
        # TODO: possibly allow Complex too
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "controller PE")
        # Link to the other side of the conversion
        p.add(cs(cst.OTHER_SIDE), "input PE", "Conversion", "output PE")
        # Make sure the two sides are not the same
        p.add(_bpp('constraint.Equality')(False), "input PE", "output PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "output PE", "output simple PE")
        # Link to ER
        p.add(cb.peToER(), "output simple PE", "output simple ER")
        p.add(_bpp('constraint.Equality')(True), "input simple ER",
              "output simple ER")
        # Make sure the input/output is a Complex
        p.add(tp(_bpimpl('Complex')().getModelInterface()), "output PE")
        p.add(tp(_bpimpl('Complex')().getModelInterface()), "input PE")
        return p

    def get_gef(self):
        """Extract Gef INDRA Statements from the BioPAX model.

        This method uses a custom BioPAX Pattern
        (one that is not implemented PatternBox) to query for controlled
        BiochemicalReactions in which the same protein is in complex with
        GDP on the left hand side and in complex with GTP on the
        right hand side. This implies that the controller is a GEF for the
        GDP/GTP-bound protein.
        """
        p = self._gef_gap_base()
        s = _bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [_match_to_array(m) for m in res.toArray()]
        for r in res_array:
            controller_pe = r[p.indexOf('controller PE')]
            input_pe = r[p.indexOf('input PE')]
            input_spe = r[p.indexOf('input simple PE')]
            output_pe = r[p.indexOf('output PE')]
            output_spe = r[p.indexOf('output simple PE')]
            reaction = r[p.indexOf('Conversion')]
            control = r[p.indexOf('Control')]

            # Make sure the GEF is not a complex
            # TODO: it could be possible to extract certain complexes here, for
            # instance ones that only have a single protein
            if _is_complex(controller_pe):
                continue

            members_in = self._get_complex_members(input_pe)
            members_out = self._get_complex_members(output_pe)
            if not (members_in and members_out):
                continue
            # Make sure the outgoing complex has exactly 2 members
            # TODO: by finding matching proteins on either side, in principle
            # it would be possible to find Gef relationships in complexes
            # with more members
            if len(members_out) != 2:
                continue
            # Make sure complex starts with GDP that becomes GTP
            gdp_in = False
            for member in members_in:
                if isinstance(member, Agent) and member.name == 'GDP':
                    gdp_in = True
            gtp_out = False
            for member in members_out:
                if isinstance(member, Agent) and member.name == 'GTP':
                    gtp_out = True
            if not (gdp_in and gtp_out):
                continue
            ras_list = self._get_agents_from_entity(input_spe)
            gef_list = self._get_agents_from_entity(controller_pe)
            ev = self._get_evidence(control)
            for gef, ras in itertools.product(_listify(gef_list),
                                               _listify(ras_list)):
                st = Gef(gef, ras, evidence=ev)
                st_dec = decode_obj(st, encoding='utf-8')
                self.statements.append(st_dec)

    def get_gap(self):
        """Extract Gap INDRA Statements from the BioPAX model.

        This method uses a custom BioPAX Pattern
        (one that is not implemented PatternBox) to query for controlled
        BiochemicalReactions in which the same protein is in complex with
        GTP on the left hand side and in complex with GDP on the
        right hand side. This implies that the controller is a GAP for the
        GDP/GTP-bound protein.
        """
        p = self._gef_gap_base()
        s = _bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [_match_to_array(m) for m in res.toArray()]
        for r in res_array:
            controller_pe = r[p.indexOf('controller PE')]
            input_pe = r[p.indexOf('input PE')]
            input_spe = r[p.indexOf('input simple PE')]
            output_pe = r[p.indexOf('output PE')]
            output_spe = r[p.indexOf('output simple PE')]
            reaction = r[p.indexOf('Conversion')]
            control = r[p.indexOf('Control')]

            # Make sure the GAP is not a complex
            # TODO: it could be possible to extract certain complexes here, for
            # instance ones that only have a single protein
            if _is_complex(controller_pe):
                continue

            members_in = self._get_complex_members(input_pe)
            members_out = self._get_complex_members(output_pe)
            if not (members_in and members_out):
                continue
            # Make sure the outgoing complex has exactly 2 members
            # TODO: by finding matching proteins on either side, in principle
            # it would be possible to find Gap relationships in complexes
            # with more members
            if len(members_out) != 2:
                continue
            # Make sure complex starts with GDP that becomes GTP
            gtp_in = False
            for member in members_in:
                if isinstance(member, Agent) and member.name == 'GTP':
                    gtp_in = True
            gdp_out = False
            for member in members_out:
                if isinstance(member, Agent) and member.name == 'GDP':
                    gdp_out = True
            if not (gtp_in and gdp_out):
                continue
            ras_list = self._get_agents_from_entity(input_spe)
            gap_list = self._get_agents_from_entity(controller_pe)
            ev = self._get_evidence(control)
            for gap, ras in itertools.product(_listify(gap_list),
                                               _listify(ras_list)):
                st = Gap(gap, ras, evidence=ev)
                st_dec = decode_obj(st, encoding='utf-8')
                self.statements.append(st_dec)

    @staticmethod
    def _get_complex_members(cplx):
        # Get the members of a complex. This is returned as a list
        # of lists since complexes can contain other complexes. The
        # list of lists solution allows us to preserve this.
        member_pes = cplx.getComponent().toArray()

        # Make a dict of member URIs and their
        # corresponding stoichiometries
        member_stos = {s.getPhysicalEntity().getUri():
                        s.getStoichiometricCoefficient() for
                        s in cplx.getComponentStoichiometry().toArray()}

        # Some complexes do not have any members explicitly listed
        if not member_pes:
            member_pes = cplx.getMemberPhysicalEntity().toArray()
            if not member_pes:
                logger.debug('Complex "%s" has no members.' %
                            cplx.getDisplayName())
                return None
        members = []
        for m in member_pes:
            if _is_complex(m):
                ms = BiopaxProcessor._get_complex_members(m)
                if ms is None:
                    return None
                members.extend(ms)
            else:
                ma = BiopaxProcessor._get_agents_from_entity(m)
                try:
                    sto = member_stos[m.getUri()]
                    sto_int = int(sto)
                except KeyError:
                    # No stoichiometry information - assume it is 1
                    members.append(ma)
                    sto_int = 1
                for i in range(sto_int):
                    members.append(ma)
        return members

    @staticmethod
    def _get_entity_mods(bpe):
        """Get all the modifications of an entity in INDRA format"""
        if _is_entity(bpe):
            features = bpe.getFeature().toArray()
        else:
            features = bpe.getEntityFeature().toArray()
        mods = []
        for feature in features:
            if not _is_modification(feature):
                continue
            mc = BiopaxProcessor._extract_mod_from_feature(feature)
            if mc is not None:
                mods.append(mc)
        return mods

    def _get_generic_modification(self, mod_class):
        """Get all modification reactions given a Modification class."""
        mod_type = modclass_to_modtype[mod_class]
        if issubclass(mod_class, RemoveModification):
            mod_gain_const = mcct.LOSS
            mod_type = modtype_to_inverse[mod_type]
        else:
            mod_gain_const = mcct.GAIN
        mod_filter = mod_type[:5]
        # Start with a generic modification pattern
        p = BiopaxProcessor._construct_modification_pattern()
        p.add(mcc(mod_gain_const, mod_filter),
                  "input simple PE", "output simple PE")
        s = _bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [_match_to_array(m) for m in res.toArray()]
        stmts = []
        for r in res_array:
            controller_pe = r[p.indexOf('controller PE')]
            input_pe = r[p.indexOf('input PE')]
            input_spe = r[p.indexOf('input simple PE')]
            output_spe = r[p.indexOf('output simple PE')]
            reaction = r[p.indexOf('Conversion')]
            control = r[p.indexOf('Control')]

            if not _is_catalysis(control):
                continue
            cat_dir = control.getCatalysisDirection()
            if cat_dir is not None and cat_dir.name() != 'LEFT_TO_RIGHT':
                logger.debug('Unexpected catalysis direction: %s.' % \
                    control.getCatalysisDirection())
                continue

            enzs = BiopaxProcessor._get_primary_controller(controller_pe)
            if not enzs:
                continue
            '''
            if _is_complex(input_pe):
                sub_members_in = self._get_complex_members(input_pe)
                sub_members_out = self._get_complex_members(output_pe)
                # TODO: It is possible to find which member of the complex is
                # actually modified. That member will be the substrate and
                # all other members of the complex will be bound to it.
                logger.info('Cannot handle complex substrates.')
                continue
            '''
            subs = BiopaxProcessor._get_agents_from_entity(input_spe,
                                                           expand_pe=False)
 
            ev = self._get_evidence(control)
            for enz, sub in itertools.product(_listify(enzs), _listify(subs)):
                # Get the modifications
                mod_in = \
                    BiopaxProcessor._get_entity_mods(input_spe)
                mod_out = \
                    BiopaxProcessor._get_entity_mods(output_spe)

                sub.mods = _get_mod_intersection(mod_in, mod_out)

                if issubclass(mod_class, AddModification):
                    gained_mods = _get_mod_difference(mod_out, mod_in)
                else:
                    gained_mods = _get_mod_difference(mod_in, mod_out)

                for mod in gained_mods:
                    # Is it guaranteed that these are all modifications
                    # of the type we are extracting?
                    if mod.mod_type not in (mod_type,
                                            modtype_to_inverse[mod_type]):
                        continue
                    stmt = mod_class(enz, sub, mod.residue, mod.position,
                                     evidence=ev)
                    stmts.append(decode_obj(stmt, encoding='utf-8'))
        return stmts

    @staticmethod
    def _get_primary_controller(controller_pe):
        # If it's not a complex, just return the corresponding agent
        if not _is_complex(controller_pe):
            enzs = BiopaxProcessor._get_agents_from_entity(controller_pe)
            return enzs

        # Identifying the "real" enzyme in a complex may not always be
        # possible.
        # One heuristic here could be to find the member which is
        # active and if it is the only active member then
        # set this as the enzyme to which all other members of the
        # complex are bound.
        # Get complex members
        members = BiopaxProcessor._get_complex_members(controller_pe)
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

    @staticmethod
    def _construct_modification_pattern():
        """Construct the BioPAX pattern to extract modification reactions."""
        # The following constraints were pieced together based on the
        # following two higher level constrains: pb.controlsStateChange(),
        # pb.controlsPhosphorylation().
        p = _bpp('Pattern')(_bpimpl('PhysicalEntity')().getModelInterface(),
                            'controller PE')
        # Getting the control itself
        p.add(cb.peToControl(), "controller PE", "Control")
        # Link the control to the conversion that it controls
        p.add(cb.controlToConv(), "Control", "Conversion")
        # The controller shouldn't be a participant of the conversion
        p.add(_bpp('constraint.NOT')(cb.participant()),
              "Conversion", "controller PE")
        # Get the input participant of the conversion
        p.add(pt(rt.INPUT, True), "Control", "Conversion", "input PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "input PE", "input simple PE")
        # Link to ER
        p.add(cb.peToER(), "input simple PE", "input simple ER")
        # Make sure the participant is a protein
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "input simple PE")
        # Link to the other side of the conversion
        p.add(cs(cst.OTHER_SIDE), "input PE", "Conversion", "output PE")
        # Make sure the two sides are not the same
        p.add(_bpp('constraint.Equality')(False), "input PE", "output PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "output PE", "output simple PE")
        # Link to ER
        p.add(cb.peToER(), "output simple PE", "output simple ER")
        p.add(_bpp('constraint.Equality')(True), "input simple ER",
              "output simple ER")
        # Make sure the output is a Protein
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "output simple PE")
        p.add(_bpp('constraint.NOT')(cb.linkToSpecific()),
              "input PE", "output simple PE")
        p.add(_bpp('constraint.NOT')(cb.linkToSpecific()),
              "output PE", "input simple PE")
        return p

    @staticmethod
    def _get_agent_from_entity(bpe):
        name = BiopaxProcessor._get_element_name(bpe)
        db_refs = BiopaxProcessor._get_db_refs(bpe)
        if _is_protein(bpe):
            mcs = BiopaxProcessor._get_entity_mods(bpe)
        else:
            mcs = []
        agent = Agent(name, db_refs=db_refs, mods=mcs)
        return agent

    @staticmethod
    def _get_agents_from_entity(bpe, expand_pe=True, expand_er=True):
        # If the entity has members (like a protein family),
        # we iterate over them
        if expand_pe:
            members = bpe.getMemberPhysicalEntity().toArray()
            if members:
                agents = []
                for m in members:
                    member_agents = BiopaxProcessor._get_agents_from_entity(m)
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
                members = er.getMemberEntityReference().toArray()
                if members:
                    agents = []
                    for m in members:
                        agent = BiopaxProcessor._get_agent_from_entity(m)
                        # For entity references, we remove context
                        agent.mods = []
                        agents.append(agent)
                    return agents
        # If it is a single entity, we get its name and database
        # references
        agent = BiopaxProcessor._get_agent_from_entity(bpe)
        return agent

    @staticmethod
    def _extract_mod_from_feature(mf):
        """Extract the type of modification and the position from
        a ModificationFeature object in the INDRA format."""
        # ModificationFeature / SequenceModificationVocabulary
        mf_type = mf.getModificationType()
        if mf_type is None:
            return None
        mf_type_terms = mf_type.getTerm().toArray()
        known_mf_type = None
        for t in mf_type_terms:
            if t.startswith('MOD_RES '):
                t = t[8:]
            mf_type_indra = _mftype_dict.get(t)
            if mf_type_indra is not None:
                known_mf_type = mf_type_indra
                break
        if not known_mf_type:
            logger.debug('Skipping modification with unknown terms: %s' %
                        ', '.join(mf_type_terms))
            return None

        mod_type, residue = known_mf_type

        # getFeatureLocation returns SequenceLocation, which is the
        # generic parent class of SequenceSite and SequenceInterval.
        # Here we need to cast to SequenceSite in order to get to
        # the sequence position.
        mf_pos = mf.getFeatureLocation()
        if mf_pos is not None:
            # If it is not a SequenceSite we can't handle it
            if not mf_pos.modelInterface.getName() == \
                'org.biopax.paxtools.model.level3.SequenceSite':
                mod_pos = None
            else:
                mf_site = cast(_bp('SequenceSite'), mf_pos)
                mf_pos_status = mf_site.getPositionStatus()
                if mf_pos_status is None:
                    mod_pos = None
                elif mf_pos_status and mf_pos_status.toString() != 'EQUAL':
                    logger.debug('Modification site position is %s' %
                                mf_pos_status.toString())
                else:
                    mod_pos = mf_site.getSequencePosition()
                    mod_pos = '%s' % mod_pos
        else:
            mod_pos = None
        mc = ModCondition(mod_type, residue, mod_pos, True)
        return mc

    @staticmethod
    def _get_evidence(bpe):
        citations = BiopaxProcessor._get_citations(bpe)
        source_id = bpe.getUri()
        if not citations:
            citations = [None]
        epi = {'direct': True}
        sources = bpe.getDataSource().toArray()
        annotations = {}
        if sources:
            if len(sources) > 1:
                logger.warning('More than one data source for %s' % bpe.uri)
            if sources[0].uri:
                entry = sources[0].uri.split('/')[-1]
                annotations['source_sub_id'] = entry
        ev = [Evidence(source_api='biopax', pmid=cit,
                       source_id=source_id, epistemics=epi,
                       annotations=annotations)
              for cit in citations]
        return ev

    @staticmethod
    def _get_citations(bpe):
        xrefs = bpe.getXref().toArray()
        refs = []
        for xr in xrefs:
            db_name = xr.getDb()
            if db_name is not None and db_name.upper() == 'PUBMED':
                refs.append(xr.getId())
        # TODO: handle non-pubmed evidence
        # TODO: do we need to look at bpe.getEvidence()
        return refs

    @staticmethod
    def _get_db_refs(bpe):
        db_refs = {}
        if _is_protein(bpe) or _is_rna(bpe):
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            # Handle missing HGNC/UP ids
            if hgnc_id and not uniprot_id:
                uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
            elif uniprot_id and not hgnc_id:
                if uniprot_client.is_human(uniprot_id):
                    hgnc_name = uniprot_client.get_gene_name(uniprot_id, False)
                    if hgnc_name:
                        hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
            # If we have both an HGNC ID and a Uniprot ID, override the 
            # Uniprot ID with the one associated with the HGNC ID
            elif uniprot_id and hgnc_id:
                hgnc_up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if hgnc_up_id != uniprot_id:
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
    def _get_element_name(bpe):
        def get_name(bpe):
            # FIXME Deal with case when HGNC entry is not name
            # Deal with case when multiple Uniprot IDs marked as
            # primary
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            if hgnc_id is not None:
                name = BiopaxProcessor._get_hgnc_name(hgnc_id)
                if name is None:
                    name = bpe.getDisplayName()
            elif uniprot_id is not None:
                name = uniprot_client.get_gene_name(uniprot_id)
                if name is None:
                    name = bpe.getDisplayName()
            else:
                name = bpe.getDisplayName()
            return name

        if _is_protein(bpe) or _is_rna(bpe):
            name = get_name(bpe)
        elif _is_small_molecule(bpe):
            name = bpe.getDisplayName()
        elif _is_physical_entity(bpe):
            name = bpe.getDisplayName()
        else:
            logger.debug('Unhandled entity type %s' %
                        bpe.getModelInterface().getName())
            name = bpe.getDisplayName()

        return name

    @staticmethod
    def _get_uniprot_id(bpe):
        # There is often more than one UniProt ID reported.
        # This usually corresponds to the primary accession ID and one or more
        # secondary accession IDs (these IDs are from deprecated entries that
        # have been merged into the primary.
        def map_to_up_primary(ids):
            primary_ids = []
            for up_id in ids:
                if not uniprot_client.is_secondary(up_id):
                    primary_ids.append(up_id)
                    continue
                primary_id = uniprot_client.get_primary_id(up_id)
                primary_ids.append(primary_id)
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
                                'choosing the first: %s'
                                 % ','.join(primary_ids))
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
        uri = bp_entref.getUri()
        # First try to match the URI itself to see if it is a UniProt
        # reference.
        m = re.match('http://identifiers.org/uniprot/([A-Z0-9]+)', uri)
        if m:
            uniprot_id = m.groups()[0]
            primary_id = map_to_up_primary([uniprot_id])
            return primary_id
        # If the URI is not a UniProt reference then we look through xrefs
        xrefs = bp_entref.getXref().toArray()
        uniprot_refs = [x for x in xrefs if
                        (x.getDb() is not None and
                         x.getDb().lower() in ('uniprot knowledgebase',
                                               'uniprotkb'))]
        if not uniprot_refs:
            return None
        uniprot_ids = [r.getId() for r in uniprot_refs]
        primary_id = map_to_up_primary(uniprot_ids)
        return primary_id

    @staticmethod
    def _get_chebi_id(bpe):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.getXref().toArray()
        chebi_ids = []
        for xr in xrefs:
            dbname = xr.getDb()
            dbid = xr.getId()
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
                                 (dbid, bpe.getDisplayName()))
        if not chebi_ids:
            return None
        elif len(chebi_ids) == 1:
            return chebi_ids[0]
        else:
            return chebi_ids

    @staticmethod
    def _get_rna_grounding(bpe):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return {}
        xrefs = bp_entref.getXref().toArray()
        rna_grounding = {}
        for xr in xrefs:
            dbname = xr.getDb()
            dbid = xr.getId()
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
    def _get_chemical_grounding(bpe):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return {}
        xrefs = bp_entref.getXref().toArray()
        chemical_grounding = {}
        for xr in xrefs:
            dbname = xr.getDb()
            dbid = xr.getId()
            if dbname is None:
                continue
            dbname = dbname.upper()
            if dbname == 'PUBCHEM-COMPOUND':
                chemical_grounding['PUBCHEM'] = 'PUBCHEM:%s' % dbid
            elif dbname == 'MESH':
                chemical_grounding['MESH'] = dbid
            elif dbname == 'DRUGBANK':
                chemical_grounding['DRUGBANK'] = dbid
            elif dbname == 'HMDB':
                chemical_grounding['HMDB'] = dbid
        return chemical_grounding

    @staticmethod
    def _get_hgnc_id(bpe):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.getXref().toArray()
        # Check for HGNC IDs
        hgnc_ids = [x.getId() for x in xrefs if
                    (x.getDb() is not None and x.getDb().lower() == 'hgnc')]
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
            hgnc_syms = [x.getId() for x in xrefs
                         if (x.getDb() is not None and
                             x.getDb().lower() == 'hgnc symbol')]
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
    def _get_hgnc_name(hgnc_id):
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        return hgnc_name

    @staticmethod
    def _get_entref(bpe):
        """Returns the entity reference of an entity if it exists or
        return the entity reference that was passed in as argument."""
        if not _is_reference(bpe):
            try:
                er = bpe.getEntityReference()
            except AttributeError:
                return None
            return er
        else:
            return bpe

    def get_coverage(self):
        uids = set()
        objs = self.model.getObjects()
        for obj in objs.toArray():
            if isinstance(obj, _bpimpl('Catalysis')) or \
                isinstance(obj, _bpimpl('TemplateReactionRegulation')):
                uids.add(obj.getUri())
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

# Functions for accessing frequently used java classes with shortened path
def _bp(path):
    prefix = 'org.biopax.paxtools.model.level3'
    classname = prefix + '.' + path
    return _autoclass_robust(classname)


def _bpp(path):
    prefix = 'org.biopax.paxtools.pattern'
    classname = prefix + '.' + path
    return _autoclass_robust(classname)


def _bpimpl(path):
    prefix = 'org.biopax.paxtools.impl.level3'
    postfix = 'Impl'
    classname = prefix + '.' + path + postfix
    return _autoclass_robust(classname)


def _autoclass_robust(path):
    try:
        cl = autoclass(path)
    except JavaException:
        logger.error('Could not instantiate ' + path)
        return None
    return cl


def _cast_biopax_element(bpe):
    """ Casts a generic BioPAXElement object into a specific type.
    This is useful when a search only returns generic elements. """
    return cast(bpe.getModelInterface().getName(), bpe)

def _match_to_array(m):
    """ Returns an array consisting of the elements obtained from a pattern
    search cast into their appropriate classes. """
    return [_cast_biopax_element(m.get(i)) for i in range(m.varSize())]

def _is_complex(pe):
    """Return True if the physical entity is a complex"""
    val = isinstance(pe, _bp('Complex')) or \
            isinstance(pe, _bpimpl('Complex'))
    return val

def _is_protein(pe):
    """Return True if the element is a protein"""
    val = isinstance(pe, _bp('Protein')) or \
            isinstance(pe, _bpimpl('Protein')) or \
            isinstance(pe, _bp('ProteinReference')) or \
            isinstance(pe, _bpimpl('ProteinReference'))
    return val

def _is_rna(pe):
    """Return True if the element is an RNA"""
    val = isinstance(pe, _bp('Rna')) or isinstance(pe, _bpimpl('Rna'))
    return val

def _is_small_molecule(pe):
    """Return True if the element is a small molecule"""
    val = isinstance(pe, _bp('SmallMolecule')) or \
            isinstance(pe, _bpimpl('SmallMolecule')) or \
            isinstance(pe, _bp('SmallMoleculeReference')) or \
            isinstance(pe, _bpimpl('SmallMoleculeReference'))
    return val

def _is_physical_entity(pe):
    """Return True if the element is a physical entity"""
    val = isinstance(pe, _bp('PhysicalEntity')) or \
            isinstance(pe, _bpimpl('PhysicalEntity'))
    return val

def _is_modification(feature):
    return (_is_modification_or_activity(feature) == 'modification')

def _is_activity(feature):
    return (_is_modification_or_activity(feature) == 'activity')

def _is_modification_or_activity(feature):
    """Return True if the feature is a modification"""
    if not (isinstance(feature, _bp('ModificationFeature')) or \
            isinstance(feature, _bpimpl('ModificationFeature'))):
        return None
    mf_type = feature.getModificationType()
    if mf_type is None:
        return None
    mf_type_terms = mf_type.getTerm().toArray()
    for term in mf_type_terms:
        if term in ('residue modification, active',
                    'residue modification, inactive',
                    'active', 'inactive'):
            return 'activity'
    return 'modification'

def _is_reference(bpe):
    """Return True if the element is an entity reference."""
    if isinstance(bpe, _bp('ProteinReference')) or \
        isinstance(bpe, _bpimpl('ProteinReference')) or \
        isinstance(bpe, _bp('SmallMoleculeReference')) or \
        isinstance(bpe, _bpimpl('SmallMoleculeReference')) or \
        isinstance(bpe, _bp('RnaReference')) or \
        isinstance(bpe, _bpimpl('RnaReference')) or \
        isinstance(bpe, _bp('EntityReference')) or \
        isinstance(bpe, _bpimpl('EntityReference')):
        return True
    else:
        return False

def _is_entity(bpe):
    """Return True if the element is a physical entity."""
    if isinstance(bpe, _bp('Protein')) or \
        isinstance(bpe, _bpimpl('Protein')) or \
        isinstance(bpe, _bp('SmallMolecule')) or \
        isinstance(bpe, _bpimpl('SmallMolecule')) or \
        isinstance(bpe, _bp('Complex')) or \
        isinstance(bpe, _bpimpl('Complex')) or \
        isinstance(bpe, _bp('Rna')) or \
        isinstance(bpe, _bpimpl('Rna')) or \
        isinstance(bpe, _bp('RnaRegion')) or \
        isinstance(bpe, _bpimpl('RnaRegion')) or \
        isinstance(bpe, _bp('DnaRegion')) or \
        isinstance(bpe, _bpimpl('DnaRegion')) or \
        isinstance(bpe, _bp('PhysicalEntity')) or \
        isinstance(bpe, _bpimpl('PhysicalEntity')):
        return True
    else:
        return False

def _is_catalysis(bpe):
    """Return True if the element is Catalysis."""
    if isinstance(bpe, _bp('Catalysis')) or \
        isinstance(bpe, _bpimpl('Catalysis')):
        return True
    else:
        return False

def _has_members(bpe):
    if _is_reference(bpe):
        members =  bpe.getMemberEntityReference().toArray()
    elif _is_entity(bpe):
        members =  bpe.getMemberPhysicalEntity().toArray()
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

def _get_mod_intersection(mods1, mods2):
    shared_mods = []
    for m1 in mods1:
        found = False
        for m2 in mods2:
            if m1.matches(m2):
                found = True
                break
        if found:
            shared_mods.append(m1)
    return shared_mods

def _get_mod_difference(mods1, mods2):
    difference_mods = []
    for m1 in mods1:
        found = False
        for m2 in mods2:
            if m1.matches(m2):
                found = True
                break
        if not found:
            difference_mods.append(m1)
    return difference_mods


# Some BioPAX Pattern classes as shorthand
pb = _bpp('PatternBox')
cb = _bpp('constraint.ConBox')
rt = _bpp('util.RelType')
tp = _bpp('constraint.Type')
cs = _bpp('constraint.ConversionSide')
cst = _bpp('constraint.ConversionSide$Type')
pt = _bpp('constraint.Participant')
mcc = _bpp('constraint.ModificationChangeConstraint')
mcct = _bpp('constraint.ModificationChangeConstraint$Type')
