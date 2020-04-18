import networkx
from indra.databases import hgnc_client
from indra.databases import uniprot_client
from indra.databases import chebi_client
from indra.databases import mesh_client


class IndraOntology(networkx.MultiDiGraph):
    def __init__(self):
        super().__init__()
        self.add_hgnc_nodes()
        self.add_uniprot_nodes()
        self.add_hgnc_uniprot_links()

    def add_hgnc_nodes(self):
        nodes = [(label('HGNC', hid),
                  {'name': hname, 'ns': 'HGNC', 'id': hid})
                 for (hid, hname) in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        nodes = [(label('UP', uid), {'name': uname, 'ns': 'UP', 'id': uid})
                 for (uid, uname)
                 in uniprot_client.um.uniprot_gene_name.items()]
        self.add_nodes_from(nodes)

    def add_hgnc_uniprot_links(self):
        edges = []
        for hid, uid in hgnc_client.uniprot_ids.items():
            uids = uid.split(', ')
            for uid in uids:
                edges.append((label('HGNC', hid), label('UP', uid),
                              {'type': 'xref', 'source': 'hgnc'}))
        self.add_edges_from(edges)

        edges = [(label('UP', uid), label('HGNC', hid),
                  {'type': 'xref', 'source': 'hgnc'})
                 for uid, hid in uniprot_client.um.uniprot_hgnc.items()]
        self.add_edges_from(edges)


def label(ns, id):
    return '%s:%s' % (ns, id)
