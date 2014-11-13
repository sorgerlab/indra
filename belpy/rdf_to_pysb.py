import re
import sys
import keyword
import urllib

import rdflib
from rdflib import URIRef, Namespace
from rdflib.namespace import RDF

from pysb import *
from pysb import SelfExporter, InvalidComponentNameError

SelfExporter.do_export = False

"""
For example graphs showing how BEL is represented in RDF, see:

http://wiki.openbel.org/display/BEL2RDF/BEL

Documentation for rdflib can be found at

https://rdflib.readthedocs.org

Types of uncertainty
--------------------

- Uncertainty about initial conditions
    - Concentrations of protein and mRNA species
    - Activity states of signaling molecules (e.g., how much AKT is "active"
      in basal conditions?)
- Uncertainty about kinetic rate parameters
- Uncertainty about structural aspects of interactions (binding sites,
  modification sites)
    - Which modifications on kinases are the activating ones?

Todo
----
- Kin -> Kin rules
- Phosphatase --> Kin rules
- Get all complexes and make binding rules
- Get amino acid substitutions ("sub" terms in protein abundances)

Abundance types
---------------

ProteinAbundance
ModifiedProteinAbundance
AbundanceActivity

NOTE:
ModifiedProteinAbundances will be a subset of ProteinAbundances, since all
ModifiedProteinAbundances have a protein abundance as a child (via the
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
    # Replace any spaces or hyphens with underscores
    term = term.replace(' ', '_')
    term = term.replace('-', '_')
    return term


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
    gene_names = [name_from_uri(r[0]) for r in res_prots]

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
        protein = name_from_uri(mod[0])
        mod_type = term_from_uri(mod[1])
        pos = term_from_uri(mod[2])
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
        protein = name_from_uri(act[0])
        act_type = term_from_uri(act[1])
        agents[protein]['Activity'].add(act_type)

    # ASSEMBLY -----
    # Now we assemble the PySB model.
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
            state_value = ['active', 'inactive']
            site_list.append(act)
            state_dict[state_key] = state_value

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

def get_statements(g, model):
    # Query for all statements where a kinase directlyIncreases modified
    # form of substrate. Ignore kinase activity of complexes for now and
    # include only the kinase activities of ProteinAbundances.
    q_phospho = prefixes + """
        SELECT ?kinaseName ?substrateName ?mod ?pos ?subject ?object
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

    for stmt in res_phospho:
        kin_name = name_from_uri(stmt[0])
        sub_name = name_from_uri(stmt[1])
        mod = term_from_uri(stmt[2])
        mod_pos = term_from_uri(stmt[3])
        # For the rule names: unfortunately, due to what looks like a bug in
        # the BEL to RDF conversion, the statements themselves are stringified
        # as, e.g.,
        # http://www.openbel.org/bel/kin_p_HGNC_KDR_http://www.openbel.org/vocabulary/DirectlyIncreases_p_HGNC_KDR_pmod_P_Y_996
        # where the subject and object are separated by a URI-prefixed relationship
        # term. This screws up the term_from_uri function which strips the
        # URI off. As a result I've manually reconstituted valid names here.
        subj = term_from_uri(stmt[4])
        obj = term_from_uri(stmt[5])
        rule_name = name_from_uri('%s_directlyIncreases_%s' % (subj, obj))
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

        rule = Rule(rule_name,
                    kin_mono(Kinase='active') + sub_mono(**{site_name: 'u'}) >>
                    kin_mono(Kinase='active') + sub_mono(**{site_name: 'p'}),
                    kf_phospho)
        model.add_component(rule)

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
    model = get_monomers(g)
    get_statements(g, model)
