import os
import re
import sys
import ipdb
import warnings
from pysb import *
from pysb import ReactionPattern, ComplexPattern, ComponentDuplicateNameError

# The CLASSPATH environmental variable needs to be set
# to point to the paxtools jar file.
# This can be obtained from
# http://sourceforge.net/projects/biopax/files/paxtools/
import jnius_config
jnius_config.add_options('-Xmx4g')
from jnius import autoclass, JavaException


site_states = {
    'phospho': ['u', 'p'],
    'activity': ['inactive', 'active'],
    'sumo': ['n', 'y'],
    'methyl': ['n', 'y'],
    'glyco': ['n', 'y'],
    'ub': ['n', 'y'],
    'acetyl': ['n', 'y'],
    'hydroxyl': ['n', 'y']
    }

residues = {
    'serine': 'S',
    'threonine': 'T',
    'tyrosine': 'Y',
    'cysteine': 'C',
    'lysine': 'K',
    'proline': 'P',
    'arginine': 'R',
    'asparagine': 'N'
    }

def iscomplex(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.ComplexImpl"))
def isprotein(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.ProteinImpl"))
def issmallmolecule(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.SmallMoleculeImpl"))
def ismodfeature(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.ModificationFeatureImpl"))
def iscatalysis(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.CatalysisImpl"))
def iscontrol(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.ControlImpl"))
def isreaction(e):
    return isinstance(e, autoclass("org.biopax.paxtools.impl.level3.BiochemicalReactionImpl"))


def get_valid_name(name):
    name = name.replace(' ', '_')
    name = name.replace('-', '_')
    name = name.replace('/', '_')
    name = name.replace('(', '')
    name = name.replace(')', '')
    name = name.replace('+', '_plus')
    name = name.replace('-', '_minus')
    if not re.match(r'[a-z]\Z', name[0], re.IGNORECASE):
        name = 'p' + name
    return name

def get_residue_name(residue):
    try:
        residue_name = residues[residue]
    except KeyError:
        warnings.warn('Unknown residue %s' % residue)
        return ''
    return residue_name

def get_modification_site(f):
    site_acc = autoclass("org.biopax.paxtools.controller.PathAccessor") \
        ("ModificationFeature/featureLocation:SequenceSite/sequencePosition")
    site = site_acc.getValueFromBean(f).toArray()
    if site:
        if len(site) > 1:
            warnings.warn('More than one site')
        site = site[0]
        site_name = '%d' % site
        return site_name
    else:
        return ''

def get_modification_pattern(feature_array):
    pattern = {}
    for f in feature_array:
        try:
            mod_type, mod_loc, mod_state = get_modification_feature(f)
        except TypeError:
            continue
        mod_name = mod_type
        if mod_loc:
            mod_name += '_' + mod_loc
        pattern[mod_name] = mod_state
    return pattern

def get_modification_feature(f):
    term_acc = autoclass("org.biopax.paxtools.controller.PathAccessor") \
        ("ModificationFeature/modificationType/term")
    term = term_acc.getValueFromBean(f).toArray()

    if not term:
        warnings.warn("Empty terms field")
        return None
    elif len(term) > 1:
        warnings.warn("Terms has %d elements" % len(term))
    term = term[0]

    # Activity
    if term == 'residue modification, active':
        return ('activity', '', 'active')
    if term == 'residue modification, inactive':
        return ('activity', '', 'inactive')

    if term.find('phospho') != -1:
        if term == 'phosphorylated residue':
            return ('phospho', '', 'p')
        else:
            residue = term.split('-')[-1]
            residue_name = get_residue_name(residue)
            site_name = get_modification_site(f)
            return ('phospho', residue_name + site_name, 'p')
    if term.find('acetylated'):
        residue = term.split('-')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('acetyl', residue_name + site_name, 'y')
    if term.find('hydroxylated'):
        residue = term.split(' ')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('hydroxyl', residue_name + site_name, 'y')
    if term.find('sumoylated'):
        residue = term.split(' ')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('sumo', residue_name + site_name, 'y')
    if term.find('methylated'):
        residue = term.split('-')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('methyl', residue_name + site_name, 'y')
    if term.find('glycosyl'):
        residue = term.split('-')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('glycosyl', residue_name + site_name, 'y')
    if term.find('ubiquitin'):
        if term == 'ubiquitination':
            return ('ub', '', 'y')
        if term.find('ubiquitinylated'):
            residue = term.split(' ')[-1]
            residue_name = get_residue_name(residue)
            site_name = get_modification_site(f)
            return ('ub', residue_name + site_name, 'y')
    if term.find('farnesyl'):
        residue = term.split('-')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('farnesyl', residue_name + site_name, 'y')

    if term.find('palmitoyl'):
        residue = term.split('-')[-1]
        residue_name = get_residue_name(residue)
        site_name = get_modification_site(f)
        return ('palm', residue_name + site_name, 'y')

    warnings.warn('Unhandled modification %s' % term)
    return None

def create_site(monomer, site, states=None):
    if site not in monomer.sites:
        monomer.sites.append(site)
    if states is not None:
        monomer.site_states.setdefault(site, [])
        for state in states:
            if state not in monomer.site_states[site]:
                monomer.site_states[site].append(state)

def get_create_monomer(model, biopax_element):
    if isprotein(biopax_element) or issmallmolecule(biopax_element):
        name = biopax_element.getDisplayName()
        name = get_valid_name(name)
        monomer = model.monomers.get(name)
        if monomer is None:
            monomer  = Monomer(name)
            model.add_component(monomer)
        features = biopax_element.getFeature().toArray()
        for f in features:
            try:
                mod_type, mod_loc, _ = get_modification_feature(f)
            except TypeError:
                continue
            try:
                mod_states = site_states[mod_type]
            except KeyError:
                warnings.warn('Modification %s not found in site states' % mod_type)
                continue
            if mod_loc:
               mod_type += '_' + mod_loc
            create_site(monomer, mod_type, mod_states)
    # Call recursively for each element of a complex
    elif iscomplex(biopax_element):
        for e in biopax_element.getComponent().toArray():
            create_monomer(model, e)
    else:
        warnings.warn("Create_monomer for this type of biopax element is not implemented")
        return None
    return monomer

def get_reaction_players(reaction):
    left = reaction.getLeft().toArray()
    right = reaction.getRight().toArray()
    catalyzers = []
    controllers = []
    directions = []
    if reaction.getControlledOf().size() > 0:
        catalyses = [x for x in reaction.getControlledOf().toArray() if iscatalysis(x)]
        controls = [x for x in reaction.getControlledOf().toArray() if iscontrol(x)]
        for c in catalyses:
            catalyzers += c.getController().toArray()
        for c in controls:
            controllers += c.getController().toArray()
        for c in catalyses:
            direction = c.getCatalysisDirection()
            if direction is not None:
                directions.append(direction.name())
            else:
                directions.append('')
    return (left, right, catalyzers, controllers, directions)

def get_pattern_string(monomer, pattern):
    s = monomer.name
    for site, state in pattern.iteritems():
        s += '_%s_%s' % (site, state)
    return s

def check_equal_pattern(m1, m2):
    if m1.monomer == m2.monomer:
        if m1.site_conditions == m2.site_conditions:
            return True
    return False

def assemble(model, reaction):
    left, right, catalyzers, controllers, directions = get_reaction_players(reaction)
    if len(left) != 1:
        warnings.warn('Cannot handle %d elements on left hand side.' % len(left))
        return
    if len(right) != 1:
        warnings.warn('Cannot handle %d elements on right hand side.' % len(right))
        return
    if iscomplex(left[0]) or iscomplex(right[0]):
        warnings.warn('Cannot handle complexes.')
        return
    left_monomer = get_create_monomer(model, left[0])
    if left_monomer is None:
        return
    left_pattern = get_modification_pattern(left[0].getFeature().toArray())
    right_monomer = get_create_monomer(model, right[0])
    if right_monomer is None:
        return
    right_pattern = get_modification_pattern(right[0].getFeature().toArray())

    # Break it up for one rule per catalyzer
    # TODO: handle controllers which have types like INHIBITION, etc.
    for c, d in zip(catalyzers, directions):
        if iscomplex(c):
            warnings.warn('Cannot handle complexes.')
            continue
        c_monomer = get_create_monomer(model, c)
        if c_monomer is None:
            continue
        c_pattern = get_modification_pattern(c.getFeature().toArray())

        if d == 'RIGHT_TO_LEFT':
            tmp = left_monomer
            left_monomer = right_monomer
            right_monomer = tmp
            tmp = left_pattern
            left_pattern = right_pattern
            right_pattern = tmp

        rule_name = '%s_%s_%s' % (get_pattern_string(left_monomer, left_pattern),
                                  get_pattern_string(c_monomer, c_pattern),
                                  get_pattern_string(right_monomer, right_pattern))
        kf_default = model.parameters['kf_default']
        if check_equal_pattern(left_monomer(left_pattern), right_monomer(right_pattern)):
            warnings.warn("Both sides of the reaction are identical.")
            continue

        try:
            r = Rule(rule_name,
                    left_monomer(left_pattern) + c_monomer(c_pattern) >>
                    right_monomer(right_pattern) + c_monomer(c_pattern),
                    kf_default)
            model.add_component(r)
        except ComponentDuplicateNameError:
            warnings.warn("Rule %s already in model." % rule_name)


class BiopaxProcessor(object):
    def __init__(self, biopax_model):
        self.biopax_model = biopax_model
        self.reactions = []
    def get_reactions(self):
        for r in biopax_model.getObjects().toArray():
            if isreaction(r):
                self.reactions.append(r)
    def make_model(self):
        model = Model()
        kf_default = Parameter('kf_default', 1.0)
        model.add_component(kf_default)
        for r in self.reactions:
            assemble(model, r)
        return model

def biopax_to_pysb(biopax_model):
    bp = BiopaxProcessor(biopax_model)
    bp.get_reactions()
    model = bp.make_model()
    return model


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: python biopax_to_pysb.py file.owl"
        sys.exit()

    data_file = sys.argv[1]

    io_class = autoclass('org.biopax.paxtools.io.SimpleIOHandler')
    io = io_class(autoclass('org.biopax.paxtools.model.BioPAXLevel').L3)

    try:
        fileIS = autoclass('java.io.FileInputStream')(data_file)
    except JavaException:
        print 'Could not open data file %s' % data_file
        sys.exit(0)
    try:
        biopax_model = io.convertFromOWL(fileIS)
    except JavaException:
        print 'Could not convert data file %s to BioPax model' % data_file
        sys.exit(0)

    fileIS.close()

    model = biopax_to_pysb(biopax_model)
