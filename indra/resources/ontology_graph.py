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

    def add_hgnc_nodes(self):
        nodes = [('HGNC:%s' % hid, {'name': hname}) for (hid, hname)
                 in hgnc_client.hgnc_names.items()]
        self.add_nodes_from(nodes)

    def add_uniprot_nodes(self):
        nodes = [('UP:%s' % uid, {'name': uname}) for (uid, uname)
                 in uniprot_client.um.uniprot_gene_name.items()]
        self.add_nodes_from(nodes)
