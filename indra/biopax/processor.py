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
from indra.databases import hgnc_client, uniprot_client
from indra.statements import *
from indra.biopax import pathway_commons_client as pcc
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
        """Extract INDRA Complex statements from the model."""
        for obj in self.model.getObjects().toArray():
            bpe = _cast_biopax_element(obj)
            if not _is_complex(bpe):
                continue
            ev = self._get_evidence(bpe)

            members = self._get_complex_members(bpe)
            if members is not None:
                if len(members) > 10:
                    logger.info('Skipping complex with more than 10 members.')
                    continue
                complexes = _get_combinations(members)
                for c in complexes:
                    self.statements.append(decode_obj(Complex(c, ev),
                                                      encoding='utf-8'))

    def get_modifications(self):
        """Extract INDRA Modification statements from model."""
        for modclass, modtype in modclass_to_modtype.items():
            # TODO: we could possibly try to also extract generic
            # modifications here
            if modtype == 'modification':
                continue
            stmts = self._get_generic_modification(modclass)
            self.statements += stmts

    def get_activity_modification(self):
        """Extract INDRA ActiveForm statements from the model."""
        mcc = _bpp('constraint.ModificationChangeConstraint')
        mcct = _bpp('constraint.ModificationChangeConstraint$Type')
        mod_filter = 'residue modification, active'
        for is_active in [True, False]:
            p = self._construct_modification_pattern()
            if is_active:
                rel = mcct.GAIN
            else:
                rel = mcct.LOSS
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

    def get_regulate_amounts(self):
        """Extract INDRA RegulateAmount statements from the model."""
        pb = _bpp('PatternBox')

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
                logger.info('Complex "%s" has no members.' %
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
        feats = [f for f in bpe.getFeature().toArray() if _is_modification(f)]
        mods = []
        for f in feats:
            mc = BiopaxProcessor._extract_mod_from_feature(f)
            if mc is not None:
                mods.append(mc)
        return mods

    def _get_generic_modification(self, mod_class):
        """Get all modification reactions given a Modification class."""
        mcc = _bpp('constraint.ModificationChangeConstraint')
        mcct = _bpp('constraint.ModificationChangeConstraint$Type')
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
                logger.info('Unexpected catalysis direction: %s.' % \
                    control.getCatalysisDirection())
                continue
            if _is_complex(controller_pe):
                # Identifying the "real" enzyme in a complex may not always be
                # possible.
                # One heuristic here could be to find the member which is
                # active and if it is the only active member then
                # set this as the enzyme to which all other members of the
                # complex are bound.
                # Get complex members
                members = self._get_complex_members(controller_pe)
                if members is None:
                    continue
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
                            logger.info(msg)
                            continue
                else:
                    msg = 'Cannot handle complex enzymes with ' + \
                            'multiple protein members.'
                    logger.info(msg)
                    continue
            else:
                enzs = BiopaxProcessor._get_agents_from_entity(controller_pe)
            if _is_complex(input_pe):
                # It is possible to find which member of the complex is
                # actually modified. That member will be the substrate and
                # all other members of the complex will be bound to it.
                logger.info('Cannot handle complex substrates.')
                continue
            ev = self._get_evidence(control)


            subs = BiopaxProcessor._get_agents_from_entity(input_spe,
                                                           expand_pe=False)
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
    def _construct_modification_pattern():
        '''
        Constructs the BioPAX pattern to extract modification reactions
        '''
        pb = _bpp('PatternBox')
        cb = _bpp('constraint.ConBox')
        flop = _bpp('constraint.Field$Operation')
        rt = _bpp('util.RelType')
        tp = _bpp('constraint.Type')
        cs = _bpp('constraint.ConversionSide')
        cst = _bpp('constraint.ConversionSide$Type')
        pt = _bpp('constraint.Participant')

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
        p.add(_bpp('constraint.Equality')(True), "input simple ER", "output simple ER")
        # Make sure the output is a Protein
        p.add(tp(_bpimpl('Protein')().getModelInterface()), "output simple PE")
        p.add(_bpp('constraint.NOT')(cb.linkToSpecific()),
              "input PE", "output simple PE")
        p.add(_bpp('constraint.NOT')(cb.linkToSpecific()),
              "output PE", "input simple PE")
        return p

    @staticmethod
    def _get_agent_from_entity(bpe):
        if bpe.getMemberPhysicalEntity().toArray():
            raise ValueError('Not a single entity.')
        name = BiopaxProcessor._get_element_name(bpe)
        db_refs = BiopaxProcessor._get_db_refs(bpe)
        mods = BiopaxProcessor._get_entity_mods(bpe, get_activity=False)
        mcs = [ModCondition(m[0], m[1], m[2], True) for m in mods]
        agent = Agent(name, db_refs=db_refs, mods=mcs)
        
        

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
        mcs = BiopaxProcessor._get_entity_mods(bpe)

        if expand_er:
            er = BiopaxProcessor._get_entref(bpe)
            if er is not None:
                members = er.getMemberEntityReference().toArray()
                if members:
                    agents = []
                    for m in members:
                        name = BiopaxProcessor._get_element_name(m)
                        db_refs = BiopaxProcessor._get_db_refs(m)
                        agents.append(Agent(name, db_refs=db_refs, mods=mcs))
                    return agents
        # If it is a single entity, we get its name and database
        # references
        name = BiopaxProcessor._get_element_name(bpe)
        db_refs = BiopaxProcessor._get_db_refs(bpe)
        agent = Agent(name, db_refs=db_refs, mods=mcs)
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
            mf_type_indra = _mftype_dict.get(t)
            if mf_type_indra is None:
                logger.info('Unknown modification type term: %s' % t)
            else:
                known_mf_type = mf_type_indra
        if not known_mf_type:
            logger.info('Ignored modification with terms: %s' %
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
                    logger.info('Modification site position is %s' %
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
        ev = [Evidence(source_api='biopax', pmid=cit,
                       source_id=source_id)
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
        elif _is_small_molecule(bpe):
            chebi_id = BiopaxProcessor._get_chebi_id(bpe)
            if chebi_id is not None:
                db_refs['CHEBI'] = chebi_id
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
            logger.info('Unhandled entity type %s' %
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
            if xr.getDb().upper() == 'CHEBI':
                chebi_ids.append(xr.getId().replace('CHEBI:', ''))
            elif xr.getDb().upper() == 'CAS':
                # Special handling of common entities
                if xr.getId() == '86-01-1':
                    chebi_ids.append('15996')
                elif xr.getId() == '24696-26-2':
                    chebi_ids.append('17761')
                else:
                    logging.info('Unknown cas id: %s' % xr.getId())
        if not chebi_ids:
            return None
        elif len(chebi_ids) == 1:
            return chebi_ids[0]
        else:
            return chebi_ids

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
    'residue modification, active': ('active', None),
    'residue modification, inactive': ('inactive', None),
    'N4-glycosyl-L-asparagine': ('glycosylation', 'N'),
    'n4glycoasn': ('glycosylation', 'N'),
    'O-glycosyl-L-threonine': ('glycosylation', 'T'),
    'S-palmitoyl-L-cysteine': ('palmitoylation', 'C'),
    'N6-acetyl-L-lysine' : ('acetylation', 'K'),
    'n6aclys': ('acetylation', 'K'),
    'naclys': ('acetylation', 'K'),
    'N-acetylated L-lysine': ('acetylation', 'K'),
    'hydroxylated proline': ('hydroxylation', 'P'),
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
                    'residue modification, inactive'):
            return 'activity'
    return 'modification'

def _is_reference(bpe):
    """Return True if the element is an entity reference."""
    if isinstance(bpe, _bp('ProteinReference')) or \
        isinstance(bpe, _bpimpl('ProteinReference')) or \
        isinstance(bpe, _bp('SmallMoleculeReference')) or \
        isinstance(bpe, _bpimpl('SmallMoleculeReference')) or \
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
