import tqdm
import pandas
from indra.statements import *
from indra.databases import hgnc_client
from indra.ontology.standardize import standardize_db_refs, \
    standardize_name_db_refs


rel_mapping = {
    'increases^expression': IncreaseAmount,
    'decreases^expression': DecreaseAmount,
    'increases^activity': Activation,
    'decreases^activity': Inhibition,
}


class CTDProcessor:
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_chemical_gene(self):
        for _, row in tqdm.tqdm(self.df.iterrows()):
            chem_name, chem_mesh_id, _, gene_name, gene_entrez_id, \
                gene_type, organism_name, organism_tax_id, txt, \
                rels, pmids = list(row)

            chem_agent = get_chemical_agent(chem_name, chem_mesh_id)
            gene_agent = get_gene_agent(gene_name, gene_entrez_id)
            stmt_types = get_statement_types(rels)
            context = get_context(organism_name, organism_tax_id)
            evs = [Evidence(pmid=pmid, text=txt, context=context)
                   for pmid in pmids.split('|')]
            for stmt_type in stmt_types:
                stmt = stmt_type(chem_agent, gene_agent, evidence=evs)
                self.statements.append(stmt)


def get_context(organism_name, organism_tax_id):
    if pandas.isna(organism_tax_id):
        return None
    tax_id = str(int(organism_tax_id))
    return RefContext(organism_name,
                      db_refs={'TAXONOMY': tax_id})


def get_statement_types(rel_str):
    rels = rel_str.split('|')
    return [rel_mapping[rel] for rel in rels if rel in rel_mapping]


def get_gene_agent(name, gene_entrez_id):
    db_refs = {'EGID': gene_entrez_id}
    hgnc_id = hgnc_client.get_hgnc_id(name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
    standard_name, db_refs = standardize_name_db_refs(db_refs)
    if standard_name:
        name = standard_name
    return Agent(name, db_refs=db_refs)


def get_chemical_agent(name, mesh_id):
    db_refs = standardize_db_refs({'MESH': mesh_id})
    return Agent(name, db_refs=db_refs)