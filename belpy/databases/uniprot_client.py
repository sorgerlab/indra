import rdflib
import urllib2

uniprot_url = 'http://www.uniprot.org/uniprot/'

rdf_prefixes = """
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX db: <http://purl.uniprot.org/database/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> """


def query_protein(protein_id):
    url = uniprot_url + protein_id + '.rdf'
    g = rdflib.Graph()
    try:
        g.parse(url)
    except urllib2.HTTPError:
        print 'Could not find protein with id %s' % protein_id
        return None
    return g


def get_hgnc_name(g):
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?res a up:Resource .
            ?res up:database db:HGNC .
            ?res rdfs:comment ?name .
            }
        """
    res = g.query(query)
    if len(res) > 0:
        hgnc_name = [r for r in res][0][0].toPython()
        return hgnc_name
    else:
        return None


def get_gene_name(g):
    # This is an alternative to get_hgnc_name and is useful when
    # HGNC name is not availabe (for instance, when the organism 
    # is not homo sapiens)
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?gene a up:Gene .
            ?gene skos:prefLabel ?name .
            }
        """
    res = g.query(query)
    if len(res) > 0:
        gene_name = [r for r in res][0][0].toPython()
        return gene_name
    else:
        return None
   

def get_sequence(g):
    query = rdf_prefixes + """
        SELECT ?seq
        WHERE {
            ?simple_seq a up:Simple_Sequence .
            ?simple_seq rdf:value ?seq .
            }
        """
    res = g.query(query)
    seq = [r for r in res][0][0].toPython()
    return seq


def get_modifications(g):
    query = rdf_prefixes + """
        SELECT ?beg_pos ?comment
        WHERE {
            ?mod_res a up:Modified_Residue_Annotation .
            ?mod_res rdfs:comment ?comment .
            ?mod_res up:range ?range .
            ?range faldo:begin ?beg .
            ?range faldo:end ?end .
            ?beg a faldo:ExactPosition .
            ?beg faldo:position ?beg_pos .
            FILTER (?beg = ?end)
            }
        """
    res = g.query(query)
    mods = []
    for r in res:
        mods.append((r[0].value, r[1].value))
    return mods

def verify_location(g, residue, location):
    """
    Verify if a given residue is at the given location
    acording to the UniProt sequence
    """
    seq = get_sequence(g)
    try:
        if seq[location] == residue:
            return True
        else:
            return False
    except IndexError:
        return False


if __name__ == '__main__':
    g = query_protein('Q02750')
    seq = get_sequence(g)
    mods = get_modifications(g)
    hgnc_name = get_hgnc_name(g)
    gene_name = get_gene_name(g)
