import ipdb
import sys

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


def print_result_generic(res):
    for r in res.toArray():
        for bpe in res.getVariables():
            print bpe.toString()
            print '================'


class BiopaxProcessor(object):
    def __init__(self, model):
        self.model = model
        self.belpy_stmts = []
        self.hgnc_cache = {}

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
        mcc = bpp('constraint.ModificationChangeConstraint')
        mcct = bpp('constraint.ModificationChangeConstraint$Type')
        # Start with a generic modification pattern
        p = self._construct_modification_pattern()
        # The modification type should contain "phospho"
        p.add(mcc(mcct.ANY, "phospho"), "input simple PE", "output simple PE")

        s = bpp('Searcher')
        res = s.searchPlain(self.model, p)
        res_array = [match_to_array(m) for m in res.toArray()]
        print '%d results found' % res.size()
        for r in res_array:
            enz_name = self._get_entity_names(r[p.indexOf('controller ER')])
            sub_name = self._get_entity_names(r[p.indexOf('changed generic ER')])
            stmt_str = ''
            citation = ''
            evidence = ''
            annotations = ''

            # Get the modification (s)
            outPE = r[p.indexOf('output PE')]
            # Do we need to look at EntityFeatures?
            outMF = [mf for mf in outPE.getFeature().toArray()
                     if isinstance(mf, bpimpl('ModificationFeature'))]
            mod_pos = []
            mod = []
            for mf in outMF:
                # ModificationFeature / SequenceModificationVocabulary
                mf_type = mf.getModificationType().getTerm().toArray()[0]
                if len(mf.getModificationType().getTerm().toArray()) > 1:
                    warnings.warn('More than one modification term')
                if mf_type == 'phosphorylated residue':
                    mod.append('Phosphorylation')
                elif mf_type == 'O-phospho-L-serine':
                    mod.append('PhosphorylationSerine')
                elif mf_type == 'O-phospho-L-threonine':
                    mod.append('PhosphorylationThreonine')
                elif mf_type == 'O-phospho-L-tyrosine':
                    mod.append('PhosphorylationTyrosine')
                else:
                    warnings.warn('Unknown phosphorylation type %s' % mf_type)
                    continue
                # getFeatureLocation returnr SequenceLocation, which is the
                # generic parent class of SequenceSite and SequenceInterval.
                # Here we need to cast to SequenceSite in order to get to
                # the sequence position.
                mf_pos = mf.getFeatureLocation()
                if mf_pos is not None:
                    mf_site = cast(bp('SequenceSite'), mf_pos)
                    if mf_site.getPositionStatus().toString() != 'EQUAL':
                        warnings.warn('Phosphorylation site position is %s' %
                                      mf_site.getPositionStatus().toString())
                    mod_pos.append(mf_site.getSequencePosition())
                else:
                    mod_pos.append(None)
                # Do we need to look at mf.getFeatureLocationType()?
                # It seems to be always None.
                mf_pos_type = mf.getFeatureLocationType()

            self.belpy_stmts.append(
                Phosphorylation(enz_name, sub_name,
                                mod, mod_pos, stmt_str,
                                citation, evidence, annotations))

    def print_statements(self):
        for i, stmt in enumerate(self.belpy_stmts):
            print "%s: %s" % (i, stmt)

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
            hgnc_name = self.hgnc_cache[hgnc_id]
        except KeyError:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
            self.hgnc_cache[hgnc_id] = hgnc_name
        return hgnc_name

if __name__ == '__main__':
    pass
