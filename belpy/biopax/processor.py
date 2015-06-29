import ipdb
import sys
import pickle

import jnius_config
jnius_config.add_options('-Xmx4g')
from jnius import autoclass, JavaException
from jnius import cast

from belpy.databases import hgnc_client
from belpy.statements import *


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


class BiopaxProcessor(object):
    def __init__(self, model):
        self.model = model
        self.belpy_stmts = []
        self._hgnc_cache = self._load_hgnc_cache()

    def get_complexes(self):
        pb = bpp('PatternBox')
        s = bpp('Searcher')
        p = pb.bindsTo()
        res = s.searchPlain(self.model, p)
        res_array = [match_to_array(m) for m in res.toArray()]
        print '%d results found' % res.size()
        for r in res_array:
            members = []
            # Extract first member
            members += self._get_entity_names(r[0])
            members += self._get_entity_names(r[4])
            # Modifications of first member
            feat_1 = r[1].getFeature().toArray()
            # Modifications of second member
            feat_2 = r[3].getFeature().toArray()
            self.belpy_stmts.append(Complex(members))

    def get_phosphorylation(self):
        stmts = self._get_generic_modification('phospho')
        for s in stmts:
            self.belpy_stmts.append(Phosphorylation(*s))

    def get_dephosphorylation(self):
        stmts = self._get_generic_modification('phospho', mod_gain=False)
        for s in stmts:
            self.belpy_stmts.append(Dephosphorylation(*s))

    def get_acetylation(self):
        stmts = self._get_generic_modification('acetyl')
        for s in stmts:
            self.belpy_stmts.append(Acetylation(*s))

    def get_glycosylation(self):
        stmts = self._get_generic_modification('glycosyl')
        for s in stmts:
            self.belpy_stmts.append(Glycosylation(*s))

    def get_palmiotylation(self):
        stmts = self._get_generic_modification('palmiotyl')
        for s in stmts:
            self.belpy_stmts.append(Glycosylation(*s))

    def get_activity_modification(self):
        mcc = bpp('constraint.ModificationChangeConstraint')
        mcct = bpp('constraint.ModificationChangeConstraint$Type')
        p = self._construct_modification_pattern()
        mod_filter = 'residue modification, active'
        p.add(mcc(mcct.GAIN, mod_filter),
              "input simple PE", "output simple PE")

        s = bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [match_to_array(m) for m in res.toArray()]
        print '%d results found' % res.size()
        for r in res_array:
            monomer_name = self._get_entity_names(r[p.indexOf('changed generic ER')])
            stmt_str = ''
            citation = self._get_citation(r[p.indexOf('Conversion')])
            evidence = ''
            annotations = ''
            outPE = r[p.indexOf('output PE')]
            activity = 'Activity'
            relationship = 'DirectlyIncreases' 
            mod, mod_pos = self._get_modification_site(outPE)
            if mod:
                stmt = ActivityModification(monomer_name, mod, mod_pos, 
                                            relationship, activity,
                                            stmt_str, citation, evidence, 
                                            annotations)
                self.belpy_stmts.append(stmt)

    def _get_modification_site(self, modPE):
        # Do we need to look at EntityFeatures?
        modMF = [mf for mf in modPE.getFeature().toArray()
                 if isinstance(mf, bpimpl('ModificationFeature'))]
        mod_pos = []
        mod = []

        for mf in modMF:
            # ModificationFeature / SequenceModificationVocabulary
            mf_type = mf.getModificationType().getTerm().toArray()[0]
            if len(mf.getModificationType().getTerm().toArray()) != 1:
                warnings.warn('Other than one modification term')
            try:
                mod_type = self._mftype_dict[mf_type]
            except KeyError:
                warnings.warn('Unknown modification type %s' % mf_type)
                continue

            mod.append(mod_type)

            # getFeatureLocation returns SequenceLocation, which is the
            # generic parent class of SequenceSite and SequenceInterval.
            # Here we need to cast to SequenceSite in order to get to
            # the sequence position.
            mf_pos = mf.getFeatureLocation()
            if mf_pos is not None:
                mf_site = cast(bp('SequenceSite'), mf_pos)
                mf_pos_status = mf_site.getPositionStatus()
                if mf_pos_status and mf_pos_status.toString() != 'EQUAL':
                    warnings.warn('Modification site position is %s' %
                                  mf_pos_status.toString())
                mod_pos.append(mf_site.getSequencePosition())
            else:
                mod_pos.append(None)
            # Do we need to look at mf.getFeatureLocationType()?
            # It seems to be always None.
            # mf_pos_type = mf.getFeatureLocationType()
        return mod, mod_pos


    def _get_generic_modification(self, mod_filter=None, mod_gain=True):
        mcc = bpp('constraint.ModificationChangeConstraint')
        mcct = bpp('constraint.ModificationChangeConstraint$Type')
        # Start with a generic modification pattern
        p = self._construct_modification_pattern()
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
        print '%d results found' % res.size()
        stmts = []
        for r in res_array:
            enz_name = self._get_entity_names(r[p.indexOf('controller ER')])
            sub_name = self._get_entity_names(r[p.indexOf('changed generic ER')])
            stmt_str = ''
            citation = self._get_citation(r[p.indexOf('Conversion')])
            evidence = ''
            annotations = ''

            # Get the modification (s)
            if mod_gain:
                modPE = r[p.indexOf('output PE')]
            else:
                modPE = r[p.indexOf('input PE')]

            # TODO: handle case when
            # r[p.indexOf('output simple PE')].getRDFId()
            # is not equal to r[p.indexOf('output PE')].getRDFId()
            mod, mod_pos = self._get_modification_site(modPE)
            stmts.append((enz_name, sub_name, mod, mod_pos,
                          stmt_str, citation, evidence, annotations))
        return stmts

    def print_statements(self):
        for i, stmt in enumerate(self.belpy_stmts):
            print "%s: %s" % (i, stmt)

    def _get_citation(self, bpe):
        evidence = bpe.getEvidence().toArray()
        refs = []
        for e in evidence:
            pub = e.getXref().toArray()
            for p in pub:
                if p.getDb() is None:
                    refs.append(p.getUrl().toArray())
                else:
                    refs.append('%s:%s' % (p.getDb(), p.getId()))
        return refs

    def _construct_modification_pattern(self):
        pb = bpp('PatternBox')
        cb = bpp('constraint.ConBox')
        flop = bpp('constraint.Field$Operation')
        rt = bpp('util.RelType')
        tp = bpp('constraint.Type')
        cs = bpp('constraint.ConversionSide')
        cst = bpp('constraint.ConversionSide$Type')
        pt = bpp('constraint.Participant')

        # The following constrainsts were pieced together based on the
        # following two higher level constrains: pb.controlsStateChange(),
        # pb.controlsPhosphorylation(). The pattern cannot be started
        # from EntityReference because it cannot be instantiated.
        # Therefore starting with ProteinReference as controller ER
        p = bpp('Pattern')(bpimpl('ProteinReference')().getModelInterface(),
                           "controller ER")
        # Getting the generic controller EntityReference
        p.add(cb.linkedER(True), "controller ER", "generic controller ER")
        # Getting the controller PhysicalEntity
        p.add(cb.erToPE(), "generic controller ER", "controller simple PE")
        # Getting to the complex controller PhysicalEntity
        p.add(cb.linkToComplex(), "controller simple PE", "controller PE")
        # Getting the control itself
        p.add(cb.peToControl(), "controller PE", "Control")
        # Link the control to the conversion that it controls
        p.add(cb.controlToConv(), "Control", "Conversion")
        # The controller shouldn't be a participant of the conversion
        p.add(bpp('constraint.NOT')(cb.participantER()),
              "Conversion", "controller ER")
        # Get the input participant of the conversion
        p.add(pt(rt.INPUT, True), "Control", "Conversion", "input PE")
        # Make sure the participant is a protein
        p.add(tp(bpimpl('Protein')().getModelInterface()), "input PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "input PE", "input simple PE")
        # Get the EntityReference for the converted entity
        p.add(cb.peToER(), "input simple PE", "changed generic ER")
        # Link to the other side of the conversion
        p.add(cs(cst.OTHER_SIDE), "input PE", "Conversion", "output PE")
        # Make sure the two sides are not the same
        p.add(bpp('constraint.Equality')(False), "input PE", "output PE")
        # Make sure the output is a Protein
        p.add(tp(bpimpl('Protein')().getModelInterface()), "output PE")
        # Get the specific PhysicalEntity
        p.add(cb.linkToSpecific(), "output PE", "output simple PE")
        # Link output to the converted EntityReference
        p.add(cb.peToER(), "output simple PE", "changed generic ER")
        # Get the specific converted EntityReference
        p.add(cb.linkedER(False), "changed generic ER", "changed ER")
        p.add(bpp('constraint.NOT')(cb.linkToSpecific()),
              "input PE", "output simple PE")
        p.add(bpp('constraint.NOT')(cb.linkToSpecific()),
              "output PE", "input simple PE")
        return p

    def _get_entity_names(self, bp_ent):
        names = []
        if isinstance(bp_ent, bp('Complex')):
            names += [self._get_entity_names(m) for
                      m in bp_ent.getComponent().toArray()]
        elif isinstance(bp_ent, bp('ProteinReference')) or \
                isinstance(bp_ent, bp('SmallMoleculeReference')) or \
                isinstance(bp_ent, bp('EntityReference')):
            hgnc_id = self._get_hgnc_id(bp_ent)
            if hgnc_id is None:
                hgnc_name = bp_ent.getDisplayName()
            else:
                hgnc_name = self._get_hgnc_name(hgnc_id)
            names += [hgnc_name]
        elif isinstance(bp_ent, bpimpl('Protein')) or \
                isinstance(bp_ent, bpimpl('SmallMolecule')) or \
                isinstance(bp_ent, bp('Protein')) or \
                isinstance(bp_ent, bp('SmallMolecule')):
            ref = bp_ent.getEntityReference()
            names += self._get_entity_names(ref)

        return names

    def _get_hgnc_id(self, bp_entref):
        xrefs = bp_entref.getXref().toArray()
        hgnc_refs = [x for x in xrefs if x.getDb() == 'HGNC']
        hgnc_id = None
        for r in hgnc_refs:
            try:
                hgnc_id = int(r.getId())
            except ValueError:
                continue
        return hgnc_id

    def _get_hgnc_name(self, hgnc_id):
        try:
            hgnc_name = self._hgnc_cache[hgnc_id]
        except KeyError:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
            self._hgnc_cache[hgnc_id] = hgnc_name
        return hgnc_name

    def _load_hgnc_cache(self):
        try:
            fh = open('hgnc_cache.pkl', 'rb')
        except IOError:
            return {}
        return pickle.load(fh)

    def _dump_hgnc_cache(self):
        fh = open('hgnc_cache.pkl', 'wb')
        pickle.dump(self._hgnc_cache, fh)
        fh.close()

    _mftype_dict = {
        'phosphorylated residue': 'Phosphorylation',
        'O-phospho-L-serine': 'PhosphorylationSerine',
        'O-phospho-L-threonine': 'PhosphorylationThreonine',
        'O-phospho-L-tyrosine': 'PhosphorylationTyrosine',
        'O4\'-phospho-L-tyrosine': 'PhosphorylationTyrosine'
        }

if __name__ == '__main__':
    pass
