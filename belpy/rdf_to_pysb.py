"""
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

Potential inconsistencies
=========================

- Site numbering
- Specific vs. generic sites ('Y' vs. 'Y1102')
- Activities ('Kinase' and 'Catalytic' for same enzyme)

Types of statements
===================

- Catalytic activity of X increases protein modification Y (generalized
  for phosphorylations, ubiquitinations, acetylations, etc.
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
- Mutations/substitutions that affect activity (e.g., RAS G12D)
- Modified protein -| or -> complex. Binding in the complex is conditional
  on the modification. Maybe also check if the modification ever occurs in
  the model?
- Protein modification decreases degradation rate

Notes on RDF representation
===========================

For example graphs showing how BEL is represented in RDF, see:

http://wiki.openbel.org/display/BEL2RDF/BEL

Documentation for rdflib can be found at

https://rdflib.readthedocs.org

Problems in the RDF representation
----------------------------------

For non-phosphorylation modifications, the RDF representation does not contain
the identity of the amino acid residue at the modification site, e.g.

catalyticActivity(proteinAbundance(HGNC:HIF1AN)) directlyIncreases proteinAbundance(HGNC:HIF1A,proteinModification(H,N,803))

will have no term for the asparagine (N) residue at position 803.

Similarly, for amino acid substitutions, there is no term at all describing
the substitution! The information about the substitution has to be parsed out
using a regular expression.

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

from belpy.statements import *

SelfExporter.do_export = False

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

def strip_statement(uri):
    uri = uri.replace(r'http://www.openbel.org/bel/', '')
    uri = uri.replace(r'http://www.openbel.org/vocabulary/', '')
    return uri

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
    PREFIX belns: <http://www.openbel.org/bel/namespace/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>"""

phospho_mods = [
    'PhosphorylationSerine',
    'PhosphorylationThreonine',
    'PhosphorylationTyrosine',
    'Phosphorylation',
]

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

class BelProcessor(object):
    def __init__(self, g):
        self.g = g
        self.belpy_stmts = []
        self.all_stmts = []
        self.converted_stmts = []

    def print_statements(self):
        for stmt in self.belpy_stmts:
            print stmt

    def get_modifications(self):
        # Query for all statements where a kinase directlyIncreases modified
        # form of substrate. Ignore kinase activity of complexes for now and
        # include only the kinase activities of ProteinAbundances.
        q_phospho = prefixes + """
            SELECT ?enzName ?actType ?substrateName ?mod ?pos ?subject ?object
                   ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject a belvoc:AbundanceActivity .
                ?subject belvoc:hasActivityType ?actType .
                ?subject belvoc:hasChild ?enzyme .
                ?enzyme a belvoc:ProteinAbundance .
                ?enzyme belvoc:hasConcept ?enzName .
                ?object a belvoc:ModifiedProteinAbundance .
                ?object belvoc:hasModificationType ?mod .
                ?object belvoc:hasChild ?substrate .
                ?substrate belvoc:hasConcept ?substrateName .
                OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
            }
        """

        # Now make the PySB for the phosphorylation
        res_phospho = self.g.query(q_phospho)

        for stmt in res_phospho:
            # Parse out the elements of the query
            enz_name = name_from_uri(stmt[0])
            act_type = name_from_uri(stmt[1])
            sub_name = name_from_uri(stmt[2])
            mod = term_from_uri(stmt[3])
            mod_pos = term_from_uri(stmt[4])
            subj = term_from_uri(stmt[5])
            obj = term_from_uri(stmt[6])
            stmt_str = strip_statement(stmt[7])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)

            #rule_name = get_rule_name(stmt[4], stmt[5], 'directlyIncreases')
            if act_type == 'Kinase' and mod in phospho_mods:
                self.belpy_stmts.append(
                        Phosphorylation(enz_name, sub_name, mod, mod_pos,
                                        subj, obj, stmt_str))
            elif act_type == 'Catalytic':
                if mod == 'Hydroxylation':
                    self.belpy_stmts.append(
                            Hydroxylation(enz_name, sub_name, mod, mod_pos,
                                        subj, obj, stmt_str))
                elif mod == 'Sumoylation':
                    self.belpy_stmts.append(
                            Sumoylation(enz_name, sub_name, mod, mod_pos,
                                        subj, obj, stmt_str))
                elif mod == 'Acetylation':
                    self.belpy_stmts.append(
                            Acetylation(enz_name, sub_name, mod, mod_pos,
                                        subj, obj, stmt_str))
                elif mod == 'Ubiquitination':
                    self.belpy_stmts.append(
                            Ubiquitination(enz_name, sub_name, mod, mod_pos,
                                        subj, obj, stmt_str))
                else:
                    print "Warning: Unknown modification type!"
                    print("Activity: %s, Mod: %s, Mod_Pos: %s" %
                          (act_type, mod, mod_pos))
            else:
                print "Warning: Unknown modification type!"
                print("Activity: %s, Mod: %s, Mod_Pos: %s" %
                      (act_type, mod, mod_pos))

            # For the rule names: unfortunately, due to what looks like a bug
            # in the BEL to RDF conversion, the statements themselves are
            # stringified as, e.g.,
            # http://www.openbel.org/bel/kin_p_HGNC_KDR_http://www.openbel.org/vocabulary/DirectlyIncreases_p_HGNC_KDR_pmod_P_Y_996
            # where the subject and object are separated by a URI-prefixed
            # relationship term. This screws up the term_from_uri function
            # which strips the URI off. As a result I've manually reconstituted
            # valid names here.
            # Get the monomer objects from the model
            #kin_mono = model.monomers[kin_name]
            #sub_mono = model.monomers[sub_name]
            # This represents just one (perhaps the simplest one)
            # interpretation of phosphorylation: pseudo first-order, in which
            # there is no binding between the kinase and substrate. Merely
            # sufficient to get some dynamics.  The form of the rule here is
            # dependent on the conversion of activity names (e.g., 'Kinase',
            # 'Phosphatase') directly from the RDF-ified BEL.  If alternative
            # PySB shorthands were developed for these activities this would
            # have to be modified.

            """
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
            """

    def get_dephosphorylations(self):
        # Query for all statements where a phosphatase directlyDecreases
        # modified form of substrate. Ignore kinase activity of complexes for
        # now and include only the kinase activities of ProteinAbundances.
        q_phospho = prefixes + """
            SELECT ?phosName ?substrateName ?mod ?pos ?subject ?object ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject belvoc:hasActivityType belvoc:Phosphatase .
                ?subject belvoc:hasChild ?phosphatase .
                ?phosphatase a belvoc:ProteinAbundance .
                ?phosphatase belvoc:hasConcept ?phosName .
                ?object a belvoc:ModifiedProteinAbundance .
                ?object belvoc:hasModificationType ?mod .
                ?object belvoc:hasChild ?substrate .
                ?substrate belvoc:hasConcept ?substrateName .
                OPTIONAL { ?object belvoc:hasModificationPosition ?pos . }
            }
        """

        # Now make the PySB for the phosphorylation
        res_phospho = self.g.query(q_phospho)

        for stmt in res_phospho:
            # Parse out the elements of the query
            phos_name = name_from_uri(stmt[0])
            sub_name = name_from_uri(stmt[1])
            mod = term_from_uri(stmt[2])
            mod_pos = term_from_uri(stmt[3])
            subj = term_from_uri(stmt[4])
            obj = term_from_uri(stmt[5])
            stmt_str = strip_statement(stmt[6])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.belpy_stmts.append(
                    Dephosphorylation(phos_name, sub_name, mod, mod_pos,
                                    subj, obj, stmt_str))

    def get_activating_mods(self):
        # Query for all statements where a kinase directlyIncreases modified
        # form of substrate. Ignore kinase activity of complexes for now and
        # include only the kinase activities of ProteinAbundances.
        q_mods = prefixes + """
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
        res_mods = self.g.query(q_mods)

        for stmt in res_mods:
            # Parse out the elements of the query
            kin_name = name_from_uri(stmt[0])
            mod = term_from_uri(stmt[1])
            mod_pos = term_from_uri(stmt[2])
            subj = term_from_uri(stmt[3])
            obj = term_from_uri(stmt[4])
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.belpy_stmts.append(
                    ActivatingModification(kin_name, mod, mod_pos, 'Kinase',
                                           subj, obj, stmt_str))

            """
            # Get the monomer objects from the model
            kin_mono = model.monomers[kin_name]

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
            """

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
        res_cmplx = self.g.query(q_cmplx)

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

    def get_ras_gefs(self):
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
        res_gef = self.g.query(q_gef)

        for stmt in res_gef:
            gef_name = name_from_uri(stmt[0])
            ras_name = name_from_uri(stmt[1])
            gef_activity = name_from_uri(stmt[2])
            subj = term_from_uri(stmt[3])
            obj = term_from_uri(stmt[4])
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.belpy_stmts.append(
                    RasGef(gef_name, gef_activity, ras_name,
                           subj, obj, stmt_str))

        """
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
        """

    def get_ras_gaps(self):
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
        res_gap = self.g.query(q_gap)

        for stmt in res_gap:
            gap_name = name_from_uri(stmt[0])
            ras_name = name_from_uri(stmt[1])
            gap_activity = name_from_uri(stmt[2])
            subj = term_from_uri(stmt[3])
            obj = term_from_uri(stmt[4])
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.belpy_stmts.append(
                    RasGap(gap_name, gap_activity, ras_name,
                           subj, obj, stmt_str))
            """
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
            """

    def get_kinase_kinase_rules(self):
        # Query for all statements where a kinase directlyIncreases modified
        # form of substrate. Ignore kinase activity of complexes for now and
        # include only the kinase activities of ProteinAbundances.
        q_phospho = prefixes + """
            SELECT ?kinaseName ?substrateName ?subject ?object ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject belvoc:hasActivityType belvoc:Kinase .
                ?subject belvoc:hasChild ?kinase .
                ?kinase a belvoc:ProteinAbundance .
                ?kinase belvoc:hasConcept ?kinaseName .
                ?object belvoc:hasActivityType belvoc:Kinase .
                ?object belvoc:hasChild ?substrate .
                ?substrate a belvoc:ProteinAbundance .
                ?substrate belvoc:hasConcept ?substrateName .
            }
        """

        res_phospho = self.g.query(q_phospho)

        # A default parameter object for phosphorylation
        for stmt in res_phospho:
            kin_name = name_from_uri(stmt[0])
            sub_name = name_from_uri(stmt[1])
            stmt_str = strip_statement(stmt[4])

            print "--------------------------------"
            print stmt[4]
            print("This statement says that:")
            print("%s kinase activity increases kinase activity of %s" %
                  (kin_name, sub_name))
            print "It doesn't specify the site."
            act_mods = []
            for bps in self.belpy_stmts:
                if type(bps) == ActivatingModification and \
                   bps.monomer_name == sub_name:
                    act_mods.append(bps)
            # If we know about an activation modification...
            if act_mods:
                print "However, I happen to know about the following"
                print "activating modifications for %s:" % sub_name
                for act_mod in act_mods:
                    print "    %s at %s" % (act_mod.mod, act_mod.mod_pos)

    def print_statement_coverage(self):
        """Display how many of the direct statements have been converted."""

        if not self.all_stmts:
            self.get_all_direct_statements()

        #print "--- All direct statements ----------"
        #print '\n'.join(self.all_stmts)
        #print
        print
        print "Total direct statements: %d" % len(self.all_stmts)
        print("Converted statements: %d" % len(self.converted_stmts))
        print
        print "--- Unconverted statements ---------"
        for stmt in self.all_stmts:
            if not stmt in self.converted_stmts:
                print stmt

    def get_all_direct_statements(self):
        """Get all directlyIncreases/Decreases statements in the corpus.
        Stores the results of the query in self.all_stmts.
        """
        q_stmts = prefixes + """
            SELECT ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                { ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases . }
                UNION
                { ?stmt belvoc:hasRelationship belvoc:DirectlyDecreases . }
            }
        """
        res_stmts = self.g.query(q_stmts)
        self.all_stmts = [strip_statement(stmt[0]) for stmt in res_stmts]

    def get_activating_subs(self):
        """
        p_HGNC_NRAS_sub_Q_61_K_DirectlyIncreases_gtp_p_HGNC_NRAS
        p_HGNC_KRAS_sub_G_12_R_DirectlyIncreases_gtp_p_PFH_RAS_Family
        p_HGNC_BRAF_sub_V_600_E_DirectlyIncreases_kin_p_HGNC_BRAF
        """
        q_mods = prefixes + """
            SELECT ?enzyme_name ?sub_label ?act_type ?subject ?object ?stmt
            WHERE {
                ?stmt a belvoc:Statement .
                ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
                ?stmt belvoc:hasSubject ?subject .
                ?stmt belvoc:hasObject ?object .
                ?subject a belvoc:ProteinAbundance .
                ?subject belvoc:hasConcept ?enzyme_name .
                ?subject belvoc:hasChild ?sub_expr .
                ?sub_expr rdfs:label ?sub_label .
                ?object a belvoc:AbundanceActivity .
                ?object belvoc:hasActivityType ?act_type .
                ?object belvoc:hasChild ?enzyme .
                ?enzyme a belvoc:ProteinAbundance .
                ?enzyme belvoc:hasConcept ?enzyme_name .
            }
        """

        # Now make the PySB for the phosphorylation
        res_mods = self.g.query(q_mods)

        for stmt in res_mods:
            # Parse out the elements of the query
            enz_name = name_from_uri(stmt[0])
            sub_expr = term_from_uri(stmt[1])
            act_type = term_from_uri(stmt[2])
            # Parse the WT and substituted residues from the node label.
            # Strangely, the RDF for substituted residue doesn't break the
            # terms of the BEL expression down into their meaning, as happens
            # for modified protein abundances. Instead, the substitution
            # just comes back as a string, e.g., "sub(V,600,E)". This code
            # parses the arguments back out using a regular expression.
            match = re.match('sub\(([A-Z]),([0-9]*),([A-Z])\)', sub_expr)
            if match:
                matches = match.groups()
                wt_residue = matches[0]
                position = matches[1]
                sub_residue = matches[2]
            else:
                print("Warning: Could not parse substitution expression %s" %
                      sub_expr)
                continue

            subj = term_from_uri(stmt[3])
            obj = term_from_uri(stmt[4])
            stmt_str = strip_statement(stmt[5])
            # Mark this as a converted statement
            self.converted_stmts.append(stmt_str)
            self.belpy_stmts.append(
                    ActivatingSubstitution(enz_name, wt_residue, position,
                                           sub_residue, act_type,
                                           subj, obj, stmt_str))

            """
            # Get the monomer objects from the model
            kin_mono = model.monomers[kin_name]

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
            """

    def get_ras_activities(self):
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
        res_ras = self.g.query(q_ras)

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
    bp = BelProcessor(g)
    #bp.get_activating_subs()
    bp.get_modifications()
    #bp.get_dephosphorylations()
    bp.get_activating_mods()
    #bp.get_ras_gefs()
    #bp.get_ras_gaps()
    #bp.print_statement_coverage()
    """
    #bp.get_ras_activities()
    #bp.print_statements()

    #model = get_monomers(g)
    #phos_stmts = get_phosphorylation_rules(g, model, rule_type='no_binding')
    #get_complexes(g, model)
    #get_kinase_kinase_rules(g, model)
    """

"""
--- Unconverted statements from RAS neighborhood ---------

-- Phosphatase activity --

phos_p_HGNC_DUSP4_DirectlyDecreases_p_HGNC_MAPK1_pmod_P_T_185

-- Kinase --> Kinase rules --

"It says here that the kinase activity of A increases the kinase activity of B.
Presumably this occurs through phosphorylation of B by A. At what site would
you expect this to occur? Here are the sites that I know about:"

kin_p_PFH_RAF_Family_DirectlyIncreases_kin_p_HGNC_MAP2K1
kin_p_HGNC_ARAF_DirectlyIncreases_kin_p_HGNC_MAP2K1
kin_p_HGNC_BRAF_DirectlyIncreases_kin_p_HGNC_MAP2K1
kin_p_HGNC_MAP3K3_DirectlyIncreases_kin_p_HGNC_MAP2K1
kin_p_HGNC_MAP3K1_DirectlyIncreases_kin_p_HGNC_MAP2K1

-- Should have been wrapped by a kinase activity! --

p_PFH_AKT_Family_DirectlyIncreases_p_HGNC_RAF1_pmod_P_S_259

-- Activity state increases binding --

gtp_p_HGNC_HRAS_DirectlyIncreases_complex_p_HGNC_BCL2_p_HGNC_HRAS

-- Activity of a complex --

kin_complex_NCH_AMP_Activated_Protein_Kinase_Complex_DirectlyIncreases_p_HGNC_RAF1_pmod_P_S_259

-- Being in a complex stimulates activity of member of complex --

complex_p_HGNC_BRAF_p_HGNC_PRKCE_p_HGNC_RPS6KB2_DirectlyIncreases_kin_p_HGNC_PRKCE

-- Binding of A to B inhibits B's ability to phosphorylate/activate C --

complex_p_HGNC_PEBP1_p_HGNC_RAF1_DirectlyDecreases_kin_p_HGNC_RAF1_DirectlyIncreases_p_HGNC_MAP2K1_pmod_P_S
complex_p_HGNC_PEBP1_p_HGNC_RAF1_DirectlyDecreases_kin_p_HGNC_RAF1_DirectlyIncreases_kin_p_HGNC_MAP2K1

-- Mutations --

p_HGNC_HRAS_sub_G_12_V_DirectlyIncreases_gtp_p_HGNC_HRAS
p_HGNC_NRAS_sub_Q_61_L_DirectlyIncreases_gtp_p_HGNC_NRAS
p_HGNC_KRAS_sub_G_13_D_DirectlyIncreases_gtp_p_HGNC_KRAS
p_HGNC_KRAS_sub_G_12_C_DirectlyIncreases_gtp_p_HGNC_KRAS
p_HGNC_NRAS_sub_Q_61_K_DirectlyIncreases_gtp_p_HGNC_NRAS
p_HGNC_KRAS_sub_G_12_R_DirectlyIncreases_gtp_p_PFH_RAS_Family
p_HGNC_BRAF_sub_V_600_E_DirectlyIncreases_kin_p_HGNC_BRAF

-- Translocation -> Activity --

tloc_p_HGNC_RASAL1_GOCCACC_GO_0005737_GOCCACC_GO_0005886_DirectlyDecreases_gtp_p_HGNC_KRAS
tloc_p_HGNC_RASA4_GOCCACC_GO_0005737_GOCCACC_GO_0005886_DirectlyDecreases_gtp_p_HGNC_KRAS

-- Ambiguous mechanism! --

(Binding is inhibitory to the kinase in this case, but how?)
p_HGNC_PEA15_DirectlyDecreases_kin_p_HGNC_MAPK1

p_HGNC_RAF1_DirectlyDecreases_kin_p_HGNC_ATM

p_HGNC_BRAF_pmod_P_S_43_DirectlyIncreases_gtp_p_PFH_RAS_Family_DirectlyIncreases_kin_p_HGNC_BRAF

"""
