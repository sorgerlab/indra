import re
import sys
import pickle
import warnings
import itertools
import collections

from indra.java_vm import autoclass, JavaException, cast

from indra.databases import hgnc_client
from indra.statements import *
from indra.biopax import pathway_commons_client as pcc

warnings.simplefilter("always")

# TODO:
# - Extract cellularLocation from each PhysicalEntity
# - Look at participantStoichiometry within BiochemicalReaction
# - Check whether to use Control or only Catalysis (Control might not
#   be direct)
# - Implement extracting modifications with Complex enzyme
# - Implement extracting modifications with Complex substrate


class BiopaxProcessor(object):
    def __init__(self, model):
        self.model = model
        self.statements = []

    def get_complexes(self, force_contains=None):
        for obj in self.model.getObjects().toArray():
            bpe = cast_biopax_element(obj)
            if not is_complex(bpe):
                continue
            citations = self._get_citations(bpe) 
            source_id = bpe.getUri()
            if not citations:
                ev = Evidence(source_api='biopax',
                               pmid=None,
                               source_id=source_id)
            else:
                ev = [Evidence(source_api='biopax',
                               pmid=cit,
                               source_id=source_id)
                      for cit in citations]

            members = self._get_complex_members(bpe)
            if members is not None:
                if len(members) > 10:
                    print 'Skipping complex with more than 10 members.'
                    continue
                complexes = get_combinations(members)
                for c in complexes:
                    self.statements.append(Complex(c, ev))

    def get_phosphorylation(self, force_contains=None):
        stmts = self._get_generic_modification('phospho',
                                               force_contains=force_contains)
        for s in stmts:
            self.statements.append(Phosphorylation(*s))

    def get_dephosphorylation(self, force_contains=None):
        stmts = self._get_generic_modification('phospho', mod_gain=False,
                                               force_contains=force_contains)
        for s in stmts:
            self.statements.append(Dephosphorylation(*s))

    def get_acetylation(self, force_contains=None):
        stmts = self._get_generic_modification('acetyl',
                                               force_contains=force_contains)
        for s in stmts:
            self.statements.append(Acetylation(*s))

    def get_glycosylation(self, force_contains=None):
        stmts = self._get_generic_modification('glycosyl',
                                               force_contains=force_contains)
        for s in stmts:
            self.statements.append(Glycosylation(*s))

    def get_palmitoylation(self, force_contains=None):
        stmts = self._get_generic_modification('palmitoyl',
                                               force_contains=force_contains)
        for s in stmts:
            self.statements.append(Palmitoylation(*s))

    def get_activity_modification(self, force_contains=None):
        mcc = bpp('constraint.ModificationChangeConstraint')
        mcct = bpp('constraint.ModificationChangeConstraint$Type')
        mod_filter = 'residue modification, active'
        for relationship in ['increases', 'decreases']:
            p = self._construct_modification_pattern()
            if relationship == 'increases':
                rel = mcct.GAIN
            else:
                rel = mcct.LOSS
            p.add(mcc(rel, mod_filter),
                  "input simple PE", "output simple PE")

            s = bpp('Searcher')
            res = s.searchPlain(self.model, p)
            res_array = [match_to_array(m) for m in res.toArray()]

            for r in res_array:
                reaction = r[p.indexOf('Conversion')]
                citations = self._get_citations(reaction)
                activity = 'Activity'
                input_spe = r[p.indexOf('input simple PE')]
                output_spe = r[p.indexOf('output simple PE')]

                # Get the modifications
                mod_in =\
                    BiopaxProcessor._get_entity_mods(input_spe)
                mod_out =\
                    BiopaxProcessor._get_entity_mods(output_spe)

                mod_shared = set(mod_in).intersection(set(mod_out))
                gained_mods = set(mod_out).difference(set(mod_in))

                # Here we get the evidence for the BiochemicalReaction
                source_id = reaction.getUri()
                citations = BiopaxProcessor._get_citations(reaction)
                if not citations:
                    ev = Evidence(source_api='biopax',
                                  pmid=None,
                                  source_id=source_id)
                else:
                    ev = [Evidence(source_api='biopax',
                                   pmid=cit,
                                   source_id=source_id)
                          for cit in citations]

                monomers = self._get_agents_from_entity(output_spe)
                for monomer in listify(monomers):
                    if force_contains is not None:
                        if momomer not in force_contains:
                            continue
                    static_mods =\
                        set(monomer.mods).difference(gained_mods)
                    monomer.mods = static_mods

                    mods = [m for m in gained_mods 
                            if m.mod_type not in ['active', 'inactive']]
                    if mods:
                        stmt = ActivityModification(monomer, mods,
                                                relationship, activity,
                                                evidence=ev)
                        self.statements.append(stmt)

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
                warnings.warn('Complex "%s" has no members.' %
                              cplx.getDisplayName())
                return None
        members = []
        for m in member_pes:
            if is_complex(m):
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
    def _get_entity_mods(bpe, get_activity=True):
        """Get all the modifications of an entity in INDRA format"""
        feats = [f for f in bpe.getFeature().toArray() if is_modification(f)]
        mods = []
        for f in feats:
            mc = BiopaxProcessor._extract_mod_from_feature(f)
            if mc is not None:
                if not get_activity and mc.mod_type in ['active', 'inactive']:
                    # Skip activity as a modification state for now
                    continue
                mods.append(mc)
        return mods

    def _get_generic_modification(self, mod_filter=None, mod_gain=True, 
                                  force_contains=None):
        '''
        Get all modification reactions given a filter
        '''
        mcc = bpp('constraint.ModificationChangeConstraint')
        mcct = bpp('constraint.ModificationChangeConstraint$Type')
        # Start with a generic modification pattern
        p = BiopaxProcessor._construct_modification_pattern()
        # The modification type should contain the filter string
        if mod_filter is not None:
            if mod_gain:
                mod_gain_const = mcct.GAIN
            else:
                mod_gain_const = mcct.LOSS
            p.add(mcc(mod_gain_const, mod_filter),
                      "input simple PE", "output simple PE")

        s = bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [match_to_array(m) for m in res.toArray()]
        stmts = []
        for r in res_array:
            controller_pe = r[p.indexOf('controller PE')]
            input_pe = r[p.indexOf('input PE')]
            input_spe = r[p.indexOf('input simple PE')]
            output_spe = r[p.indexOf('output simple PE')]
            reaction = r[p.indexOf('Conversion')]
            control = r[p.indexOf('Control')]

            if not is_catalysis(control):
                continue
            cat_dir = control.getCatalysisDirection()
            if cat_dir is not None and cat_dir.name() != 'LEFT_TO_RIGHT':
                warnings.warn('Unexpected catalysis direction: %s.' % control.getCatalysisDirection())
                continue
            if is_complex(controller_pe):
                # Identifying the "real" enzyme in a complex may not always be
                # possible.
                # One heuristic here could be to find the member which is
                # active and if it is the only active member then
                # set this as the enzyme to which all other members of the
                # complex are bound.
                warnings.warn('Cannot handle complex enzymes.')
                continue
            if is_complex(input_pe):
                # It is possible to find which member of the complex is 
                # actually modified. That member will be the substrate and 
                # all other members of the complex will be bound to it.
                warnings.warn('Cannot handle complex substrates.')
                continue
            # TODO: should this be the citation for the control?
            # Sometimes there is an xref within Catalysis which refers to 
            # a pubmed article in a bp:PublicationXref tag.
            citations = BiopaxProcessor._get_citations(control)
            source_id = control.getUri()
            if not citations:
                ev = Evidence(source_api='biopax',
                               pmid=None,
                               source_id=source_id)
            else:
                ev = [Evidence(source_api='biopax',
                               pmid=cit,
                               source_id=source_id)
                      for cit in citations]

            enzs = BiopaxProcessor._get_agents_from_entity(controller_pe)
            subs = BiopaxProcessor._get_agents_from_entity(input_spe,
                                                           expand_pe=False)
            for enz, sub in itertools.product(listify(enzs), listify(subs)):
                # If neither the required enzyme nor the substrate is 
                # present then skip
                if force_contains is not None:
                    if (enz.name not in force_contains) and \
                        (sub.name not in force_contains):
                        continue

                # Get the modifications
                mod_in =\
                    BiopaxProcessor._get_entity_mods(input_spe)
                mod_out =\
                    BiopaxProcessor._get_entity_mods(output_spe)

                mod_shared = set(mod_in).intersection(set(mod_out))

                sub.mods = mod_shared

                if mod_gain:
                    gained_mods = set(mod_out).difference(set(mod_in))
                else:
                    gained_mods = set(mod_in).difference(set(mod_out))

                for m in gained_mods:
                    if m.mod_type  in ['active', 'inactive']:
                        # Skip activity as a modification state
                        continue
                    stmt = (enz, sub, m.residue, m.position, ev)
                    stmts.append(stmt)
        return stmts

    @staticmethod
    def _get_citations(bpe):
        xrefs = bpe.getXref().toArray()
        refs = [x.getId() for x in xrefs if x.getDb() == 'PUBMED']
        # TODO: handle non-pubmed evidence
        return refs

    @staticmethod
    def _get_evidence(bpe):
        ev = bpe.getEvidence().toArray()
        print ev
        for e in ev:
            xrefs =  e.getXref().toArray()
            # There are also evidence codes that we could extract.
            # ev_codes = e.getEvidenceCode().toArray()
        return xrefs

    @staticmethod
    def _construct_modification_pattern():
        '''
        Constructs the BioPAX pattern to extract modification reactions
        '''
        pb = bpp('PatternBox')
        cb = bpp('constraint.ConBox')
        flop = bpp('constraint.Field$Operation')
        rt = bpp('util.RelType')
        tp = bpp('constraint.Type')
        cs = bpp('constraint.ConversionSide')
        cst = bpp('constraint.ConversionSide$Type')
        pt = bpp('constraint.Participant')

        # The following constraints were pieced together based on the
        # following two higher level constrains: pb.controlsStateChange(),
        # pb.controlsPhosphorylation().
        p = bpp('Pattern')(bpimpl('PhysicalEntity')().getModelInterface(),
                           'controller PE')
        # Getting the control itself
        p.add(cb.peToControl(), "controller PE", "Control")
        # Link the control to the conversion that it controls
        p.add(cb.controlToConv(), "Control", "Conversion")
        # The controller shouldn't be a participant of the conversion
        p.add(bpp('constraint.NOT')(cb.participant()),
              "Conversion", "controller PE")
        # Get the input participant of the conversion
        p.add(pt(rt.INPUT, True), "Control", "Conversion", "input PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "input PE", "input simple PE")
        # Link to ER
        p.add(cb.peToER(), "input simple PE", "input simple ER")
        # Make sure the participant is a protein
        p.add(tp(bpimpl('Protein')().getModelInterface()), "input simple PE")
        # Link to the other side of the conversion
        p.add(cs(cst.OTHER_SIDE), "input PE", "Conversion", "output PE")
        # Make sure the two sides are not the same
        p.add(bpp('constraint.Equality')(False), "input PE", "output PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "output PE", "output simple PE")
        # Link to ER
        p.add(cb.peToER(), "output simple PE", "output simple ER")
        p.add(bpp('constraint.Equality')(True), "input simple ER", "output simple ER")
        # Make sure the output is a Protein
        p.add(tp(bpimpl('Protein')().getModelInterface()), "output simple PE")
        p.add(bpp('constraint.NOT')(cb.linkToSpecific()),
              "input PE", "output simple PE")
        p.add(bpp('constraint.NOT')(cb.linkToSpecific()),
              "output PE", "input simple PE")
        return p

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
        mods = BiopaxProcessor._get_entity_mods(bpe)

        if expand_er:
            er = BiopaxProcessor._get_entref(bpe)
            if er is not None:
                members = er.getMemberEntityReference().toArray()
                if members:
                    agents = []
                    for m in members:
                        name = BiopaxProcessor._get_element_name(m)
                        db_refs = BiopaxProcessor._get_db_refs(m)
                        agents.append(Agent(name, db_refs=db_refs, mods=mods))
                    return agents
        # If it is a single entity, we get its name and database
        # references
        name = BiopaxProcessor._get_element_name(bpe)
        db_refs = BiopaxProcessor._get_db_refs(bpe)
        agent = Agent(name, db_refs=db_refs, mods=mods)
        return agent

    @staticmethod
    def _extract_mod_from_feature(mf):
        """Extract the type of modification and the position from
        a ModificationFeature object in the INDRA format."""
        # ModificationFeature / SequenceModificationVocabulary
        mf_type = mf.getModificationType()
        if mf_type is None:
            return None
        if len(mf_type.getTerm().toArray()) != 1:
            warnings.warn('Other than one modification term')
        mf_type = mf_type.getTerm().toArray()[0]
        try:
            mod_type, residue = BiopaxProcessor._mftype_dict[mf_type]
        except KeyError:
            warnings.warn('Ignored modification type %s' % mf_type)
            return None

        # getFeatureLocation returns SequenceLocation, which is the
        # generic parent class of SequenceSite and SequenceInterval.
        # Here we need to cast to SequenceSite in order to get to
        # the sequence position.
        mf_pos = mf.getFeatureLocation()
        if mf_pos is not None:
            mf_site = cast(bp('SequenceSite'), mf_pos)
            mf_pos_status = mf_site.getPositionStatus()
            if mf_pos_status is None:
                mod_pos = None
            elif mf_pos_status and mf_pos_status.toString() != 'EQUAL':
                warnings.warn('Modification site position is %s' %
                              mf_pos_status.toString())
            else:
                mod_pos = mf_site.getSequencePosition()
        else:
            mod_pos = None
        mc = ModCondition(mod_type, residue, mod_pos)
        return mc

    @staticmethod
    def _get_db_refs(bpe):
        db_refs = {}
        if is_protein(bpe):
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            uniprot_id = BiopaxProcessor._get_uniprot_id(bpe)
            if hgnc_id is not None:
                db_refs['HGNC'] = hgnc_id
            if uniprot_id is not None:
                db_refs['UP'] = uniprot_id
        elif is_small_molecule(bpe):
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
    def _get_element_name(bpe):
        if is_protein(bpe):
            hgnc_id = BiopaxProcessor._get_hgnc_id(bpe)
            if hgnc_id is not None:
                name = BiopaxProcessor._get_hgnc_name(hgnc_id)
                if name is None:
                    name = bpe.getDisplayName()
            else:
                name = bpe.getDisplayName()
        elif is_small_molecule(bpe):
            name = bpe.getDisplayName()
        elif is_physical_entity(bpe):
            name = bpe.getDisplayName()
        else:
            warnings.warn('Unhandled entity type %s' %
                bpe.getModelInterface().getName())
            import ipdb; ipdb.set_trace()
            name = bpe.getDisplayName()

        # Canonicalize name
        name = re.sub(r'[^\w]', '_', name)
        if re.match('[0-9]', name) is not None:
            name = 'p' + name
        return name

    @staticmethod
    def _get_uniprot_id(bpe):
        # There is often more than one UniProt ID reported.
        # This usually corresponds to the primary accession ID and one or more
        # secondary accession IDs (these IDs are from deprecated entries that
        # have been merged into the primary.
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.getXref().toArray()
        uniprot_refs = [x for x in xrefs if 
                        x.getDb() == 'UniProt Knowledgebase']
        uniprot_ids = [r.getId() for r in uniprot_refs]
        if not uniprot_ids:
            return None
        elif len(uniprot_ids) == 1:
            return uniprot_ids[0]
        else:
            return uniprot_ids

    @staticmethod
    def _get_chebi_id(bpe):
        bp_entref = BiopaxProcessor._get_entref(bpe)
        if bp_entref is None:
            return None
        xrefs = bp_entref.getXref().toArray()
        chebi_refs = [x for x in xrefs if 
                        x.getDb() == 'ChEBI']
        chebi_ids = [r.getId().replace('CHEBI:', '') for r in chebi_refs]
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
        hgnc_refs = [x for x in xrefs if x.getDb() == 'HGNC']
        hgnc_id = None
        for r in hgnc_refs:
            try:
                hgnc_id = int(r.getId())
            except ValueError:
                continue
        return hgnc_id

    @staticmethod
    def _get_hgnc_name(hgnc_id):
        hgnc_name = hgnc_client.get_hgnc_name(str(hgnc_id))
        return hgnc_name

    @staticmethod
    def _get_entref(bpe):
        """Returns the entity reference of an entity if it exists or
        return the entity reference that was passed in as argument."""
        if not is_reference(bpe):
            try:
                er = bpe.getEntityReference()
            except AttributeError:
                return None
            return er
        else:
            return bpe

    def print_statements(self):
        for i, stmt in enumerate(self.statements):
            print "%s: %s" % (i, stmt)

    def save_model(self, file_name=None):
        if file_name is None:
            print 'Missing file name'
            return
        pcc.model_to_owl(self.model, file_name)

    _mftype_dict = {
        'phosphorylated residue': ('phosphorylation', None),
        'O-phospho-L-serine': ('phosphorylation', 'S'),
        'O-phospho-L-threonine': ('phosphorylation', 'T'),
        'O-phospho-L-tyrosine': ('phosphorylation', 'Y'),
        'O4\'-phospho-L-tyrosine': ('phosphorylation', 'Y'),
        'residue modification, active': ('active', None),
        'residue modification, inactive': ('inactive', None)
        }

# Functions for accessing frequently used java classes with shortened path
def bp(path):
    prefix = 'org.biopax.paxtools.model.level3'
    classname = prefix + '.' + path
    return autoclass_robust(classname)


def bpp(path):
    prefix = 'org.biopax.paxtools.pattern'
    classname = prefix + '.' + path
    return autoclass_robust(classname)


def bpimpl(path):
    prefix = 'org.biopax.paxtools.impl.level3'
    postfix = 'Impl'
    classname = prefix + '.' + path + postfix
    return autoclass_robust(classname)


def autoclass_robust(path):
    try:
        cl = autoclass(path)
    except JavaException:
        print 'Could not instantiate ' + path
        return None
    return cl


def cast_biopax_element(bpe):
    """ Casts a generic BioPAXElement object into a specific type.
    This is useful when a search only returns generic elements. """
    return cast(bpe.getModelInterface().getName(), bpe)

def match_to_array(m):
    """ Returns an array consisting of the elements obtained from a pattern
    search cast into their appropriate classes. """
    return [cast_biopax_element(m.get(i)) for i in range(m.varSize())]

def is_complex(pe):
    """Return True if the physical entity is a complex"""
    val = isinstance(pe, bp('Complex')) or \
            isinstance(pe, bpimpl('Complex'))
    return val

def is_protein(pe):
    """Return True if the element is a protein"""
    val = isinstance(pe, bp('Protein')) or \
            isinstance(pe, bpimpl('Protein')) or \
            isinstance(pe, bp('ProteinReference')) or \
            isinstance(pe, bpimpl('ProteinReference'))
    return val

def is_small_molecule(pe):
    """Return True if the element is a small molecule"""
    val = isinstance(pe, bp('SmallMolecule')) or \
            isinstance(pe, bpimpl('SmallMolecule')) or \
            isinstance(pe, bp('SmallMoleculeReference')) or \
            isinstance(pe, bpimpl('SmallMoleculeReference'))
    return val

def is_physical_entity(pe):
    """Return True if the element is a physical entity"""
    val = isinstance(pe, bp('PhysicalEntity')) or \
           isinstance(pe, bpimpl('PhysicalEntity'))
    return val

def is_modification(fe):
    """Return True if the feature is a modification"""
    val = isinstance(fe, bp('ModificationFeature')) or \
            isinstance(fe, bpimpl('ModificationFeature'))
    return val

def is_reference(bpe):
    """Return True if the element is an entity reference."""
    if isinstance(bpe, bp('ProteinReference')) or \
        isinstance(bpe, bpimpl('ProteinReference')) or \
        isinstance(bpe, bp('SmallMoleculeReference')) or \
        isinstance(bpe, bpimpl('SmallMoleculeReference')) or \
        isinstance(bpe, bp('EntityReference')) or \
        isinstance(bpe, bpimpl('EntityReference')):
        return True
    else:
        return False

def is_entity(bpe):
    """Return True if the element is a physical entity."""
    if isinstance(bpe, bp('Protein')) or \
        isinstance(bpe, bpimpl('Protein')) or \
        isinstance(bpe, bp('SmallMolecule')) or \
        isinstance(bpe, bpimpl('SmallMolecule')) or \
        isinstance(bpe, bp('PhysicalEntity')) or \
        isinstance(bpe, bpimpl('PhysicalEntity')):
        return True
    else:
        return False

def is_catalysis(bpe):
    """Return True if the element is Catalysis."""
    if isinstance(bpe, bp('Catalysis')) or \
        isinstance(bpe, bpimpl('Catalysis')):
        return True
    else:
        return False

def has_members(bpe):
    if is_reference(bpe):
        members =  bpe.getMemberEntityReference().toArray()
    elif is_entity(bpe):
        members =  bpe.getMemberPhysicalEntity().toArray()
    else:
        return False
    if len(members) > 0:
        return True
    else:
        return False

def listify(lst):
    if not isinstance(lst, collections.Iterable):
        return [lst]
    else:
        return lst

def list_listify(lst):
    return [l if isinstance(l, collections.Iterable) else [l] for l in lst]

def get_combinations(lst):
    return itertools.product(*list_listify(lst))
