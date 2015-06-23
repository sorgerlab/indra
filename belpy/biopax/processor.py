import jnius_config
jnius_config.add_options('-Xmx4g')
from jnius import autoclass, JavaException
from jnius import cast
import ipdb

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


def modification(model, mf):
    mf = cast('org.biopax.paxtools.impl.level3.ModificationFeatureImpl', model.getByID('http://purl.org/pc2/7/' + mf))
    print mf.toString()
    p = mf.getFeatureOf().toArray()[0]
    print p.getStandardName()
    print p.getDisplayName()
    print p.getName().toArray()

def protein(model, p):
    p = cast('org.biopax.paxtools.impl.level3.ProteinImpl', model.getByID('http://purl.org/pc2/7/' + p))
    print p.getStandardName()
    print p.getDisplayName()
    print p.getName().toArray()

def owl_to_model(fname):
    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3) 
    
    try:
        fileIS = autoclass('java.io.FileInputStream')(fname)
    except JavaException:
        print 'Could not open data file %s' % fname
        sys.exit(0)
    try:
        biopax_model = io.convertFromOWL(fileIS)
    except JavaException:
        print 'Could not convert data file %s to BioPax model' % data_file
        sys.exit(0)
    
    fileIS.close()

    return biopax_model

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
    

if __name__ == '__main__':
    model = owl_to_model('BRAF.owl')
    cb = bpp('constraint.ConBox')
    pb = bpp('PatternBox')
    pc = bpp('constraint.PathConstraint')
    mcc = bpp('constraint.ModificationChangeConstraint')
    mcct = bpp('constraint.ModificationChangeConstraint$Type')
    pa = autoclass('org.biopax.paxtools.controller.PathAccessor')
    fl = bpp('constraint.Field')
    flop = bpp('constraint.Field$Operation')
    s = bpp('Searcher')
    
    # SequenceEntityReference is for everything except small molecules
    # Replace with EntityReference
    p = bpp('Pattern')(bpimpl('BiochemicalReaction')().getClass(), "r")
    p.add(pc("BiochemicalReaction/controlledOf*:Catalysis"), "r", "cat")
    p.add(pc("BiochemicalReaction/left:Protein"), "r", "lPE")
    p.add(pc("BiochemicalReaction/right:Protein"), "r", "rPE")
    p.add(pc("Catalysis/controller*:Protein"), "cat", "cPE")
    p.add(pc("Protein/entityReference:EntityReference"), "lPE", 'lER')
    p.add(pc("Protein/entityReference:EntityReference"), "rPE", 'rER')
    p.add(pc("Protein/entityReference:EntityReference"), "cPE", 'cER')
    #p.add(pc("Protein/feature*:ModificationFeature"), "lPE", "lMF")
    p.add(pc("Protein/feature*:ModificationFeature"), "rPE", "rMF")
    #p.add(pc("Protein/feature*:ModificationFeature"), "cPE", "cMF")
    #p.add(pc("ModificationFeature/modificationType*:SequenceModificationVocabulary"), "lMF", "lMFT")
    p.add(pc("ModificationFeature/modificationType*:SequenceModificationVocabulary"), "rMF", "rMFT")
    #p.add(pc("ModificationFeature/modificationType*:SequenceModificationVocabulary"), "cMF", "cMFT")
    hs = autoclass('java.util.HashSet')()
    hs.add('O-phospho-L-serine')
    hs.add('O-phospho-L-threonine')
    hs.add('O-phospho-L-tyrosine')
    #p.add(fl("SequenceModificationVocabulary/term", flop.INTERSECT, hs), 'rMFT')
    p.add(mcc(mcct.GAIN,'phospho'), 'lPE', 'rPE')
    #p.add(fl("SequenceModificationVocabulary/term", flop.NOT_INTERSECT, hs), 'lMFT')
    #p.add(fl("ModificationFeature/modificationType/term", flop.INTERSECT, hs), 'rMF')
    #p.add(fl("ModificationFeature/modificationType/term", flop.NOT_INTERSECT, hs), 'lMF')
    #p.add(pc("ModificationFeature/modificationType:SequenceModificationVocabulary"), "lmf", "lmft")
    #p.add(mcc(bpp('constraint.ModificationChangeConstraint$Type').GAIN,'phospho'))
    #p.add(fl("term", fl.Operator.INTERSECT,"O-phospho-L-serine"), "rmft")
    
    # p.add(cb.linkedER(true), "controller ER", "generic controller ER")
    # p.add(cb.erToPE(), "generic controller ER", "controller simple PE")
    # p.add(cb.linkToComplex(), "controller simple PE", "controller PE")
    # p.add(cb.peToControl(), "controller PE", "Control")
    # p.add(cb.controlToConv(), "Control", "Conversion")
    # p.add(bpp('constraint.NOT')(cb.participantER()), "Conversion", "controller ER")
    # p.add(bpp('Participant')(RelType.INPUT, true), "Control", "Conversion", "input PE")
    # p.add(cb.linkToSpecific(), "input PE", "input simple PE")
    # p.add(bpp('constraint.Type')(bpimpl('SequenceEntity')().getClass()), "input simple PE")
    # p.add(cb.peToER(), "input simple PE", "changed generic ER")
    # p.add(bpp('constraint.ConversionSide')(bpp('constraint.ConversionSide.Type').OTHER_SIDE), "input PE", "Conversion", "output PE")
    # p.add(equal(false), "input PE", "output PE")
    # p.add(cb.linkToSpecific(), "output PE", "output simple PE")
    # p.add(cb.peToER(), "output simple PE", "changed generic ER")
    # p.add(cb.linkedER(false), "changed generic ER", "changed ER")
    
    
    #p = bpp('Pattern')(bpimpl('ProteinReference')().getClass(), 'controllerPR')
    #p.add(cb.isHuman(), 'controllerPR')

    #p.add(pb.controlsStateChange(), 'controllerPR')
    #p.add(cb.erToPE(),'controllerPR', 'controllerPE')
    #p.add(cb.peToControl(),'controllerPE','control')
    #p.add(cb.controlToInter(),'control','inter')

    res = s.searchPlain(model, p)
    res_array = [match_to_array(m) for m in res.toArray()]
    print '%d results found' % res.size()
