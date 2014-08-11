import rdflib
from rdflib import URIRef
from rdflib.namespace import RDF
from rdflib import Namespace
import urllib

"""
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
consistency of this type.

Also, Histone_H3_Family have methylation at position 27 but the residue type is
not specified? (it is, but it's not in the name of the modification).

proteinAbundance(MGI:Akt1,proteinModification(P,S,473)) directlyIncreases proteinAbundance(MGI:Foxo1,proteinModification(P,T,24))

In this case, it probably should be the kinase activity of the modprot abundance.
Especially because it is "directlyIncreases", so a mechanism is being asserted.
Perhaps should detect that a phosphorylation is being referred to and require
that a kinase activity on the left hand side is being specified implicitly.

Perhaps for starters, we could assert that only activities, not abundances,
could be candidates for rules.

"""

def term_from_uri(uri):
    if uri is None:
        return None

    # Strip gene name off from URI
    gene_name = uri.rsplit('/')[-1]
    # Decode URL to handle spaces, special characters
    gene_name = urllib.unquote(gene_name)
    # Replace any spaces or hyphens with underscores
    gene_name = gene_name.replace(' ', '_')
    gene_name = gene_name.replace('-', '_')
    # If starts with a number, add an underscore FIXME
    return gene_name

BEL = Namespace("http://www.openbel.org/")

g = rdflib.Graph()
#g.parse('full_abstract1.rdf', format='nt')
g.parse('small_corpus.rdf', format='nt')
#g.parse('sample_rules.rdf', format='nt')

prefixes = """
    PREFIX belvoc: <http://www.openbel.org/vocabulary/>
    PREFIX belsc: <http://www.openbel.org/bel/>
    PREFIX belns: <http://www.openbel.org/namespace/>"""

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
}

states = {
    'PhosphorylationSerine': "['u', 'p']",
    'PhosphorylationThreonine': "['u', 'p']",
    'PhosphorylationTyrosine': "['u', 'p']",
    'Phosphorylation': "['y', 'n']",
    'Ubiquitination': "['y', 'n']",
    'Farnesylation': "['y', 'n']",
    'Hydroxylation': "['y', 'n']",
    'Acetylation': "['y', 'n']",
    'Sumoylation': "['y', 'n']",
    'Glycosylation': "['y', 'n']",
    'Methylation': "['y', 'n']",
}


def get_monomers():
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
    gene_names = [term_from_uri(r[0]) for r in res_prots]

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
        protein = term_from_uri(mod[0])
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
        protein = term_from_uri(act[0])
        act_type = term_from_uri(act[1])
        agents[protein]['Activity'].add(act_type)


    # ASSEMBLY -----
    # Now we convert to PySB monomer declarations.
    # In reality, instead of making PySB as text, would create PySB objects.
    for k, v in agents.iteritems():
        #print '%s:' % k
        #print '    Mod: %s' % v['Modification']
        #print '    Act: %s' % v['Activity']
        mstr  = "Monomer('%s', [" % k
        for mod in v['Modification']:
            if len(mod) == 2:
                mstr += "'%s%s', " % (abbrevs[mod[0]], mod[1])
            elif len(mod) == 1:
                mstr += "'%s', " % abbrevs[mod[0]]
            else:
                raise Exception("Unknown modification.")
        mstr += "], {"
        for mod in v['Modification']:
            if len(mod) == 2:
                mstr += "'%s%s': %s, " % (abbrevs[mod[0]], mod[1], states[mod[0]])
            elif len(mod) == 1:
                mstr += "'%s': %s, " % (abbrevs[mod[0]], states[mod[0]])
            else:
                raise Exception("Unknown modification.")
        mstr += "})"
        print mstr

def get_statements():
    # Query for statements
    q_stmts = prefixes + """
        SELECT ?subject ?object
        WHERE {
            ?stmt a belvoc:Statement .
            ?stmt belvoc:hasSubject ?subject .
            ?stmt belvoc:hasRelationship belvoc:DirectlyIncreases .
            ?stmt belvoc:hasObject ?object .
        }
    """

    res_stmts = g.query(q_stmts)

    prot = URIRef('http://www.openbel.org/vocabulary/ProteinAbundance')
    modprot = URIRef('http://www.openbel.org/vocabulary/ModifiedProteinAbundance')
    act = URIRef('http://www.openbel.org/vocabulary/AbundanceActivity')

    for stmt in res_stmts:
        print stmt
        sub = stmt[0]
        obj = stmt[1]
        print "-----------"
        print sub
        print obj
        for s, p, o in g.triples( (sub, RDF.type, modprot) ):
            #print "protein abundance!"
            print "found"
            print s, p, o

if __name__ == '__main__':
    get_monomers()
