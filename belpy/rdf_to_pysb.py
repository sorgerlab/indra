import re
import sys
import keyword
import urllib
import collections

import rdflib
from rdflib import URIRef, Namespace
from rdflib.namespace import RDF

from pysb import *
from pysb.core import SelfExporter, InvalidComponentNameError, \
                      ComplexPattern, ReactionPattern

SelfExporter.do_export = False

"""
For example graphs showing how BEL is represented in RDF, see:

http://wiki.openbel.org/display/BEL2RDF/BEL

Documentation for rdflib can be found at

https://rdflib.readthedocs.org

Types of uncertainty
====================

- Initial conditions
    - Concentrations of protein and mRNA species
    - Activity states of signaling molecules (e.g., how much AKT is "active"
      in basal conditions?)
- Kinetic rate parameters
- Structural aspects of interactions (which/how many binding sites,
  modification sites)
- Significance of modifications, e.g., which phospho-states are the
  activating ones?
- Significance of complexes. Is binding another protein inhibitory or
  activating?

Types of statements
===================

- Kinase -> Modified substrate (can generalize to all site-modifiers?
  Should be general for Ub, Kinase, Phosphatase, Glycos, others?
- Complex(X, Y), indicates that X and Y bind
- Modified protein -> kinase activity (can be make this general to other
  activities?)
- Kinase activity => Kinase activity. These statements are labeled as being
  "direct", i.e., mechanistic, even though the phosphorylation site on the
  substrate is unspecified.
- Activity of X -> GtpBoundActivity(Y) (i.e., X is RAS GEF)
- Activity of X -| GtpBoundActivity(Y) (i.e., X is RAS GAP)
- Substitution in X -> GtpBoundActivity(Y) (i.e., activation mutation in RAS)
- Get amino acid substitutions ("sub" terms in protein abundances)
- RasGTPases slowly hydrolyze GTP by themselves, so need to add default
  rule for each monomer with a GtpBoundActivity
- Protein families: expand out to members, or use representative?

Notes on RDF representation
===========================

Abundance types
---------------

ProteinAbundance
ModifiedProteinAbundance
AbundanceActivity

.. note::

    ModifiedProteinAbundances will be a subset of ProteinAbundances, since all
    ModifiedProteinAbundances have a ProteinAbundance as a child (via the
    relationship hasChild).

Modification types
------------------

hasModificationType (has range Activity):
Phosphorylation
PhosphorylationThreonine
PhosphorylationTyrosine
PhosphorylationSerine
Hydroxylation
Acetylation

Activity types
--------------

hasActivityType (has range Activity):
Kinase
Transcription
Phosphatase
Catalytic
GtpBound
Ubiquitination
Activity (?)

Other notes
===========

For Tsc2, we find a PhosphorylationSerine and a PhosphorylationThreonine at the
same site! So this presumably reflects an error. We need for screen for
inconsistency of this type.

Also, Histone_H3_Family have methylation at position 27 but the residue type is
not specified? (it is, but it's not in the name of the modification).

proteinAbundance(MGI:Akt1,proteinModification(P,S,473)) directlyIncreases proteinAbundance(MGI:Foxo1,proteinModification(P,T,24))

In this case, it probably should be the kinase activity of the modprot abundance.
Especially because it is "directlyIncreases", so a mechanism is being asserted.
Perhaps should detect that a phosphorylation is being referred to and require
that a kinase activity on the left hand side is being specified implicitly.

Perhaps for starters, we could assert that only activities, not abundances,
could be candidates for rules.

Redundancies
------------

What to do about the following example, in which the BEL corpus contained the
following overlapping statements? Have seen similar things where statements
regarding the family overlap those of specific genes.

Rule(u'cat_p_HGNC_SOS1_directlyIncreases_gtp_p_HGNC_KRAS',
     SOS1(Catalytic='active') + KRAS(GtpBound='inactive') >>
     SOS1(Catalytic='active') + KRAS(GtpBound='active'), kf_gef),
Rule(u'cat_p_PFH_SOS_Family_directlyIncreases_gtp_p_HGNC_NRAS',
     SOS_Family(Catalytic='active') + NRAS(GtpBound='inactive') >>
     SOS_Family(Catalytic='active') + NRAS(GtpBound='active'), kf_gef),

Protein families
----------------
When protein families appear in BEL statements, these are currently getting
converted into monomers in the resulting PySB model. Instead, these should
mapped onto representative monomers in same way. An interesting use case for
a MetaKappa inheritance-type mechanism.
"""

def name_from_uri(uri):
    """Make the URI term usable as a valid Python identifier, if possible.

    First strips of the extra URI information by calling term_from_uri,
    then checks to make sure the name is a valid Python identifier.
    Currently fixes identifiers starting with numbers by prepending with
    the letter 'p'. For other cases it raises an exception.

    This function should be called when the string that is returned is to be
    used as a PySB component name, which are required to be valid Python
    identifiers.
    """
    name = term_from_uri(uri)
    # Handle the case where the string starts with a number
    if name[0].isdigit():
        name = 'p' + name
    if re.match("[_A-Za-z][_a-zA-Z0-9]*$", name) \
            and not keyword.iskeyword(name):
        pass
    else:
        raise InvalidComponentNameError(name)

    return name

def term_from_uri(uri):
    """Basic conversion of RDF URIs to more friendly strings.

    Removes prepended URI information, and replaces spaces and hyphens with
    underscores.
    """
    if uri is None:
        return None
    # Strip gene name off from URI
    term = uri.rsplit('/')[-1]
    # Decode URL to handle spaces, special characters
    term = urllib.unquote(term)
    # Replace any spaces, hyphens, or periods with underscores
    term = term.replace(' ', '_')
    term = term.replace('-', '_')
    term = term.replace('.', '_')
    return term

def get_rule_name(subj_uri, obj_uri, relation):
    """Serializes a BEL statement to a string for use as a rule name."""
    subj = term_from_uri(subj_uri)
    obj = term_from_uri(obj_uri)
    rule_name = name_from_uri('%s_%s_%s' % (subj, relation, obj))
    return rule_name

BEL = Namespace("http://www.openbel.org/")

prefixes = """
    PREFIX belvoc: <http://www.openbel.org/vocabulary/>
    PREFIX belsc: <http://www.openbel.org/bel/>
    PREFIX belns: <http://www.openbel.org/bel/namespace/>"""

abbrevs = {
    'PhosphorylationSerine': 'S',
    'PhosphorylationThreonine': 'T',
    'PhosphorylationTyrosine': 'Y',
    'Phosphorylation': 'phospho',
    'Ubiquitination': 'ub',
    'Farnesylation': 'farnesyl',
    'Hydroxylation': 'hydroxyl',
    'Acetylation': 'acetyl',
    'Sumoylation': 'sumo',
    'Glycosylation': 'glycosyl',
    'Methylation': 'methyl',
    'Modification': 'mod',
}

active_site_names = {
    'Kinase': 'kin_site',
}

states = {
    'PhosphorylationSerine': ['u', 'p'],
    'PhosphorylationThreonine': ['u', 'p'],
    'PhosphorylationTyrosine': ['u', 'p'],
    'Phosphorylation': ['u', 'p'],
    'Ubiquitination': ['y', 'n'],
    'Farnesylation': ['y', 'n'],
    'Hydroxylation': ['y', 'n'],
    'Acetylation': ['y', 'n'],
    'Sumoylation': ['y', 'n'],
    'Glycosylation': ['y', 'n'],
    'Methylation': ['y', 'n'],
    'Modification': ['y', 'n'],
}


def get_monomers(g):
    # PROTEINS ----
    # Get the full list of proteins by querying for protein abundances
    q_prots = prefixes + """
        SELECT ?genename
        WHERE {
            ?term a belvoc:ProteinAbundance .
            ?term belvoc:hasConcept ?genename .
        }
    """

    # The result list for the proteins query
    res_prots = g.query(q_prots)

    # Parse out the gene names, we'll use these for our monomer names
    gene_names = []
    for r in res_prots:
        try:
            gene_names.append(name_from_uri(r[0]))
        except InvalidComponentNameError as e:
            print "Warning: %s" % e
            continue

    # Initialize agents dict from list of monomers.
    # For each gene name in the dict, keep a set of modifications, and a set
    # of activities--both of which are initialized to be empty:
    agents = dict.fromkeys(gene_names)
    for gene in gene_names:
        agents[gene] = {'Modification': set([]), 'Activity': set([])}

    # MODIFICATIONS ----
    # Now, query for all of the protein modifications mentioned in the corpus
    q_modprots = prefixes + """
        SELECT ?genename ?mod ?pos
        WHERE {
            ?term a belvoc:ModifiedProteinAbundance .
            ?term belvoc:hasModificationType ?mod .
            ?term belvoc:hasChild ?pterm .
            ?pterm a belvoc:ProteinAbundance .
            ?pterm belvoc:hasConcept ?genename .
            OPTIONAL { ?term belvoc:hasModificationPosition ?pos . }
        }
    """
    # The result of the modification query
    res_modprots = g.query(q_modprots)

    # For each modification returned by the query, add the modification to the
    # list for the appropriate gene
    for mod in res_modprots:
        try:
            protein = name_from_uri(mod[0])
            mod_type = term_from_uri(mod[1])
            pos = term_from_uri(mod[2])
        except InvalidComponentNameError as e:
            print "Warning: %s" % e
            continue

        if pos is not None:
            mod_type = (mod_type, pos)
        else:
            mod_type = (mod_type, )
        agents[protein]['Modification'].add(mod_type)

    # ACTIVITIES ----
    # Now query for all of the activites (e.g., kinase, transcriptional, etc.)
    q_acts = prefixes + """
        SELECT ?genename ?acttype
        WHERE {
            ?term a belvoc:AbundanceActivity .
            ?term belvoc:hasActivityType ?acttype .
            ?term belvoc:hasChild ?pterm .
            ?pterm a belvoc:ProteinAbundance .
            ?pterm belvoc:hasConcept ?genename .
        }
    """
    # The result of the activity query
    res_acts = g.query(q_acts)
    # Add the activities to the list for the appropriate gene
    for act in res_acts:
        try:
            protein = name_from_uri(act[0])
            act_type = term_from_uri(act[1])
        except InvalidComponentNameError as e:
            print "Warning: %s" % e
            continue

        agents[protein]['Activity'].add(act_type)

    # ASSEMBLY -----
    # Now that we know all the activities and modification states of the
    # proteins mentioned in the entire corpus, we assemble the PySB model.
    model = Model()
    ic_param = Parameter('default_ic', 10.)
    model.add_component(ic_param)
    for k, v in agents.iteritems():
        monomer_name = k
        site_list = []
        state_dict = {}

        # Iterate over all modifications
        for mod in v['Modification']:
            if len(mod) == 2:
                site_name = '%s%s' % (abbrevs[mod[0]], mod[1])
                state_key = site_name
                state_value = states[mod[0]]
            elif len(mod) == 1:
                site_name = '%s' % abbrevs[mod[0]]
                state_key = '%s' % abbrevs[mod[0]]
                state_value = states[mod[0]]
            else:
                raise Exception("Unknown modification.")
            site_list.append(site_name)
            state_dict[state_key] = state_value
        # Iterate over all activities
        for act in v['Activity']:
            state_key = act
            state_value = ['inactive', 'active']
            site_list.append(act)
            state_dict[state_key] = state_value
            # Add the active site for binding, if there is one
            if act in active_site_names:
                site_list.append(active_site_names[act])
        site_list.append('b') # FIXME
        m = Monomer(monomer_name, site_list, state_dict)
        model.add_component(m)
        sites = {}
        for s in m.sites:
            if s in m.site_states:
                sites[s] = m.site_states[s][0]
            else:
                sites[s] = None
        model.initial(m(sites), ic_param)
    return model

def get_phosphorylation_rules(g, model, rule_type='binding'):
    # Query for all statements where a kinase directlyIncreases modified
    # form of substrate. Ignore kinase activity of complexes for now and
    # include only the kinase activities of ProteinAbundances.
    q_phospho = prefixes + """
        SELECT ?kinaseName ?substrateName ?mod ?pos ?subject ?object ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasObject ?object .
            ?subject belvoc:hasActivityType belvoc:Kinase .
            ?subject belvoc:hasChild ?kinase .
            ?kinase a belvoc:ProteinAbundance .
            ?kinase belvoc:hasConcept ?kinaseName .
            ?object a belvoc:ModifiedProteinAbundance .
            ?object belvoc:hasModificationType ?mod .
            ?object belvoc:hasChild ?substrate .
            ?substrate belvoc:hasConcept ?substrateName .
            OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
        }
    """

    # Now make the PySB for the phosphorylation
    res_phospho = g.query(q_phospho)

    # A default parameter object for phosphorylation
    kf_phospho = Parameter('kf_phospho', 1.)
    model.add_component(kf_phospho)
    statements = []

    for stmt in res_phospho:
        kin_name = name_from_uri(stmt[0])
        sub_name = name_from_uri(stmt[1])
        mod = term_from_uri(stmt[2])
        mod_pos = term_from_uri(stmt[3])
        statements.append(stmt[6])

        # For the rule names: unfortunately, due to what looks like a bug in
        # the BEL to RDF conversion, the statements themselves are stringified
        # as, e.g.,
        # http://www.openbel.org/bel/kin_p_HGNC_KDR_http://www.openbel.org/vocabulary/DirectlyIncreases_p_HGNC_KDR_pmod_P_Y_996
        # where the subject and object are separated by a URI-prefixed relationship
        # term. This screws up the term_from_uri function which strips the
        # URI off. As a result I've manually reconstituted valid names here.
        rule_name = get_rule_name(stmt[4], stmt[5], 'directlyIncreases')
        # Get the monomer objects from the model
        kin_mono = model.monomers[kin_name]
        sub_mono = model.monomers[sub_name]
        # This represents just one (perhaps the simplest one) interpretation
        # of phosphorylation: pseudo first-order, in which there is no binding
        # between the kinase and substrate. Merely sufficient to get some dynamics.
        # The form of the rule here is dependent on the conversion of activity
        # names (e.g., 'Kinase', 'Phosphatase') directly from the RDF-ified BEL.
        # If alternative PySB shorthands were developed for these activities this
        # would have to be modified.
        if mod_pos is not None:
            site_name = '%s%s' % (abbrevs[mod], mod_pos)
        else:
            site_name = abbrevs[mod]

        if rule_type == 'no_binding':
            rule = Rule(rule_name,
                        kin_mono(Kinase='active') + sub_mono(**{site_name: 'u'}) >>
                        kin_mono(Kinase='active') + sub_mono(**{site_name: 'p'}),
                        kf_phospho)
            model.add_component(rule)
        elif rule_type == 'binding':
            rule_bind = Rule('%s_bind' % rule_name,
                        kin_mono(Kinase='active', kin_site=None) +
                        sub_mono(**{site_name: 'u'}) <>
                        kin_mono(Kinase='active', kin_site=1) %
                        sub_mono(**{site_name: ('u', 1)}),
                        kf_phospho, kf_phospho)
            rule_cat =  Rule('%s_cat' % rule_name,
                        kin_mono(Kinase='active', kin_site=1) %
                        sub_mono(**{site_name: ('u', 1)}) >>
                        kin_mono(Kinase='active', kin_site=None) +
                        sub_mono(**{site_name: 'p'}),
                        kf_phospho)
            model.add_component(rule_bind)
            model.add_component(rule_cat)
    return statements

def get_activating_mods(g, model):
    # Query for all statements where a kinase directlyIncreases modified
    # form of substrate. Ignore kinase activity of complexes for now and
    # include only the kinase activities of ProteinAbundances.
    q_phospho = prefixes + """
        SELECT ?kinaseName ?mod ?pos ?subject ?object ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasObject ?object .
            ?object belvoc:hasActivityType belvoc:Kinase .
            ?object belvoc:hasChild ?kinase .
            ?kinase a belvoc:ProteinAbundance .
            ?kinase belvoc:hasConcept ?kinaseName .
            ?subject a belvoc:ModifiedProteinAbundance .
            ?subject belvoc:hasModificationType ?mod .
            ?subject belvoc:hasChild ?kinase .
            OPTIONAL { ?subject belvoc:hasModificationPosition ?pos . }
        }
    """

    # Now make the PySB for the phosphorylation
    res_phospho = g.query(q_phospho)

    kf_activation = Parameter('kf_activation', 1e5)
    model.add_component(kf_activation)
    statements = []
    for stmt in res_phospho:
        kin_name = name_from_uri(stmt[0])
        mod = term_from_uri(stmt[1])
        mod_pos = term_from_uri(stmt[2])
        statements.append(stmt[5])
        # Get the monomer objects from the model
        kin_mono = model.monomers[kin_name]
        subj = term_from_uri(stmt[3])
        obj = term_from_uri(stmt[4])
        rule_name = name_from_uri('%s_directlyIncreases_%s' % (subj, obj))

        if mod_pos is not None:
            site_name = '%s%s' % (abbrevs[mod], mod_pos)
        else:
            site_name = abbrevs[mod]

        rule = Rule(rule_name,
                    kin_mono(**{site_name: 'p', 'Kinase': 'inactive'}) >>
                    kin_mono(**{site_name: 'p', 'Kinase': 'active'}),
                    kf_activation)
        model.add_component(rule)
    return statements

def get_complexes(g, model):
    # Query for all statements where a kinase directlyIncreases modified
    # form of substrate. Ignore kinase activity of complexes for now and
    # include only the kinase activities of ProteinAbundances.
    q_cmplx = prefixes + """
        SELECT ?term ?childName
        WHERE {
            ?term a belvoc:Term .
            ?term a belvoc:ComplexAbundance .
            ?term belvoc:hasChild ?child .
            ?child belvoc:hasConcept ?childName .
        }
    """

    # Now make the PySB for the phosphorylation
    res_cmplx = g.query(q_cmplx)

    kf_binding = Parameter('kf_binding', 1)
    model.add_component(kf_binding)

    cmplx_dict = collections.defaultdict(list)
    for stmt in res_cmplx:
        cmplx_name = term_from_uri(stmt[0])
        child_name = name_from_uri(stmt[1])
        cmplx_dict[cmplx_name].append(child_name)

    for cmplx_name, cmplx_list in cmplx_dict.iteritems():
        lhs = ReactionPattern([])
        rhs = ComplexPattern([], None)
        try:
            for monomer_name in cmplx_list:
                mono = model.monomers[monomer_name]
                mp_free = mono(b=None)
                mp_bound = mono(b=1)
                lhs = lhs + mp_free
                rhs = rhs % mp_bound
            rule_name = '%s_bind' % cmplx_name
            if not model.rules.get(rule_name):
                rule = Rule('%s_bind' % cmplx_name,
                            lhs <> rhs, kf_binding, kf_binding)
                model.add_component(rule)
        except KeyError as ke:
            print "Warning: Monomer not found, ignoring: %s" % ke

def get_gef_rules(g, model, rule_type='no_binding'):
    # First, get the statements with activities as subjects.
    q_gef = prefixes + """
        SELECT ?gefName ?rasName ?gefActivity ?subject ?object ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasObject ?object .
            ?subject a belvoc:AbundanceActivity .
            ?subject belvoc:hasActivityType ?gefActivity .
            ?subject belvoc:hasChild ?gef .
            ?gef a belvoc:ProteinAbundance .
            ?gef belvoc:hasConcept ?gefName .
            ?object a belvoc:AbundanceActivity .
            ?object belvoc:hasActivityType belvoc:GtpBound .
            ?object belvoc:hasChild ?ras .
            ?ras belvoc:hasConcept ?rasName .
        }
    """
    res_gef = g.query(q_gef)

    # A default parameter object for gef
    kf_gef = Parameter('kf_gef', 1.)
    model.add_component(kf_gef)
    statements = []
    for stmt in res_gef:
        gef_name = name_from_uri(stmt[0])
        ras_name = name_from_uri(stmt[1])
        gef_activity = name_from_uri(stmt[2])
        statements.append(stmt[5])
        rule_name = get_rule_name(stmt[3], stmt[4], 'directlyIncreases')
        # Get the monomer objects from the model
        gef_mono = model.monomers[gef_name]
        ras_mono = model.monomers[ras_name]
        if rule_type == 'no_binding':
            rule = Rule(rule_name,
                        gef_mono(**{gef_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'inactive'}) >>
                        gef_mono(**{gef_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'active'}),
                        kf_gef)
            model.add_component(rule)
    return statements

def get_gap_rules(g, model, rule_type='no_binding'):
    # First, get the statements with activities as subjects.
    q_gap = prefixes + """
        SELECT ?gapName ?rasName ?gapActivity ?subject ?object ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasObject ?object .
            ?subject a belvoc:AbundanceActivity .
            ?subject belvoc:hasActivityType ?gapActivity .
            ?subject belvoc:hasChild ?gap .
            ?gap a belvoc:ProteinAbundance .
            ?gap belvoc:hasConcept ?gapName .
            ?object a belvoc:AbundanceActivity .
            ?object belvoc:hasActivityType belvoc:GtpBound .
            ?object belvoc:hasChild ?ras .
            ?ras belvoc:hasConcept ?rasName .
        }
    """
    res_gap = g.query(q_gap)

    # A default parameter object for gap
    kf_gap = Parameter('kf_gap', 1.)
    model.add_component(kf_gap)
    statements = []
    for stmt in res_gap:
        gap_name = name_from_uri(stmt[0])
        ras_name = name_from_uri(stmt[1])
        gap_activity = name_from_uri(stmt[2])
        statements.append(stmt[5])
        rule_name = get_rule_name(stmt[3], stmt[4], 'directlyIncreases')
        # Get the monomer objects from the model
        gap_mono = model.monomers[gap_name]
        ras_mono = model.monomers[ras_name]
        if rule_type == 'no_binding':
            rule = Rule(rule_name,
                        gap_mono(**{gap_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'active'}) >>
                        gap_mono(**{gap_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'inactive'}),
                        kf_gap)
            model.add_component(rule)
    return statements

def get_all_direct_statements(g):
    q_stmts = prefixes + """
        SELECT ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            { ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases . }
            UNION
            { ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases . }
        }
    """
    res_stmts = g.query(q_stmts)
    return [stmt[0] for stmt in res_stmts]

def get_ras_rules(g, model, rule_type='no_binding'):
    # First, get the statements with activities as subjects.
    q_ras = prefixes + """
        SELECT ?stmt
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasObject ?object .
            ?subject a belvoc:AbundanceActivity .
            ?subject belvoc:hasActivityType belvoc:GtpBound .
        }
    """
    res_ras = g.query(q_ras)
    # A default parameter object for ras
    #kf_ras = Parameter('kf_ras', 1.)
    #model.add_component(kf_ras)

    for stmt in res_ras:
        print stmt[0]
    """
        ras_name = name_from_uri(stmt[0])
        ras_name = name_from_uri(stmt[1])
        ras_activity = name_from_uri(stmt[2])
        rule_name = get_rule_name(stmt[3], stmt[4], 'directlyIncreases')
        # Get the monomer objects from the model
        ras_mono = model.monomers[ras_name]
        ras_mono = model.monomers[ras_name]
        if rule_type == 'no_binding':
            rule = Rule(rule_name,
                        ras_mono(**{ras_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'active'}) >>
                        ras_mono(**{ras_activity: 'active'}) +
                        ras_mono(**{'GtpBound': 'inactive'}),
                        kf_ras)
            model.add_component(rule)
    """



if __name__ == '__main__':
    # Make sure the user passed in an RDF filename
    if len(sys.argv) < 2:
        print "Usage: python rdf_to_pysb.py file.rdf"
        sys.exit()
    # We take the RDF filename as the argument
    rdf_filename = sys.argv[1]

    # Parse the RDF
    g = rdflib.Graph()
    g.parse(rdf_filename, format='nt')
    # Build the PySB model
    all_stmts = get_all_direct_statements(g)

    model = get_monomers(g)
    phos_stmts = get_phosphorylation_rules(g, model, rule_type='no_binding')
    mod_stmts = get_activating_mods(g, model)
    #get_complexes(g, model)
    gef_stmts = get_gef_rules(g, model)
    gap_stmts = get_gap_rules(g, model)
    #get_ras_rules(g, model)
    print '\n'.join(all_stmts)
    print "Total statements: %d" % len(all_stmts)
    print("Converted statements: %d" % (len(phos_stmts) + len(mod_stmts) +
                                        len(gef_stmts) + len(gap_stmts)))
    for stmt in phos_stmts + mod_stmts + gef_stmts + gap_stmts:
        print "removing: %s" % stmt
        all_stmts.remove(stmt)
    print
    print "--- Unconverted statements ---------"
    print '\n'.join(all_stmts)

