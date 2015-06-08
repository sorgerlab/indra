import rdflib
import urllib2

uniprot_url = 'http://www.uniprot.org/uniprot/'

rdf_prefixes = """
    PREFIX up: <http://purl.uniprot.org/core/>
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

def query_sequence(g):
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

def query_modifications(g):
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
        mods.append((r[0].value,r[1].value))

    return mods

if __name__ == '__main__':
    g = query_protein('Q02750')
    seq = query_sequence(g)
    mods = query_modifications(g)
    print len([a for a in res])
