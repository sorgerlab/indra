import sys

from processor import BiopaxProcessor

from indra.java_vm import autoclass, JavaException
from indra.biopax import pathway_commons_client as pcc

def process_pc_neighborhood(gene_names, neighbor_limit=1):
    query_type = autoclass('cpath.service.GraphType').NEIGHBORHOOD
    model = pcc.run_pc_query(query_type, gene_names,
                             None, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_pc_pathsbetween(gene_names, neighbor_limit=1):
    query_type = autoclass('cpath.service.GraphType').PATHSBETWEEN
    model = pcc.run_pc_query(query_type, gene_names,
                             None, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_pc_pathsfromto(source_genes, target_genes, neighbor_limit=1):
    query_type = autoclass('cpath.service.GraphType').PATHSFROMTO
    model = pcc.run_pc_query(query_type, source_genes, 
                             target_genes, neighbor_limit)
    if model is not None:
        return process_model(model)


def process_owl(owl_filename):
    model = pcc.owl_to_model(owl_filename)
    return process_model(model)


def process_model(model):
    bproc = BiopaxProcessor(model)
    # bproc.get_complexes()
    # bproc.get_phosphorylation()
    # bproc.print_statements()
    return bproc

if __name__ == '__main__':
    pass
