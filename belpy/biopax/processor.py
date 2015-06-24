import jnius_config
jnius_config.add_options('-Xmx4g')
from jnius import autoclass, JavaException
from jnius import cast
import ipdb
import sys

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

def get_modification(bp_ent):
    bp_ent

def get_hgnc_id(bp_entref):
    xrefs = bp_entref.getXref().toArray()
    hgnc_refs = [x for x in xrefs if x.getDb() == 'HGNC']
    hgnc_id = None
    for r in hgnc_refs:
        try:
            hgnc_id = int(r.getId())
        except ValueError:
            continue
    return hgnc_id

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
            members += self.get_entity_names(r[0])
            members += self.get_entity_names(r[4])
            # Modifications of first member
            feat_1 = r[1].getFeature().toArray()
            # Modifications of second member
            feat_2 = r[3].getFeature().toArray()
            self.belpy_stmts.append(Complex(members))
            
    def get_entity_names(self, bp_ent):
        names = []
        if isinstance(bp_ent, bp('Complex')):
            names += [self.get_entity_names(m) for m in bp_ent.getComponent().toArray()]
        elif isinstance(bp_ent, bp('ProteinReference')) or \
                isinstance(bp_ent, bp('SmallMoleculeReference')) or \
                isinstance(bp_ent, bp('EntityReference')):
            hgnc_id = get_hgnc_id(bp_ent)
            if hgnc_id is None:
                hgnc_name = bp_ent.getDisplayName()
            else:
                hgnc_name = self.get_hgnc_name(hgnc_id)
            names += [hgnc_name]
        elif isinstance(bp_ent, bpimpl('Protein')) or \
                isinstance(bp_ent, bpimpl('SmallMolecule')) or \
                isinstance(bp_ent, bp('Protein')) or \
                isinstance(bp_ent, bp('SmallMolecule')):
            ref = bp_ent.getEntityReference()
            names += self.get_entity_names(ref)
        
        return names

    def get_phosphorylation(self):
        pb = bpp('PatternBox')
        s = bpp('Searcher')
        p = pb.controlsPhosphorylation()
        res = s.searchPlain(model, p)
        res_array = [match_to_array(m) for m in res.toArray()]
        print '%d results found' % res.size()

    def get_hgnc_name(self, hgnc_id):
        try:
            hgnc_name = self.hgnc_cache[hgnc_id]
        except KeyError:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
            self.hgnc_cache[hgnc_id] = hgnc_name
        return hgnc_name
    
    def print_statements(self):
        for i, stmt in enumerate(self.belpy_stmts):
            print "%s: %s" % (i, stmt)

if __name__ == '__main__':
    pass

