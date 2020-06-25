import tqdm
from indra.statements import *
from indra.databases import hgnc_client
from indra.ontology.standardize import standardize_db_refs, \
    standardize_name_db_refs


rel_mapping = {
    # Activity regulation
    'increases^activity': Activation,
    'decreases^activity': Inhibition,
    # Amount regulation
    'increases^expression': IncreaseAmount,
    'decreases^expression': DecreaseAmount,
    'increases^chemical synthesis': IncreaseAmount,
    'decreases^chemical synthesis': DecreaseAmount,
    'increases^degradation': DecreaseAmount,
    'decreases^degradation': IncreaseAmount,
    'increases^abundance': IncreaseAmount,
    'decreases^abundance': DecreaseAmount,
    'increases^stability': IncreaseAmount,
    'decreases^stability': DecreaseAmount,
    # Modification
    'increases^phosphorylation': Phosphorylation,
    'decreases^phosphorylation': Dephosphorylation,
    'increases^acetylation': Acetylation,
    'decreases^acetylation': Deacetylation,
    'increases^ubiquitination': Ubiquitination,
    'decreases^ubiquitination': Deubiquitination,
    'increases^hydroxylation': Hydroxylation,
    'decreases^hydroxylation': Dehydroxylation,
    'increases^methylation': Methylation,
    'decreases^methylation': Demethylation,
    'increases^farnesylation': Farnesylation,
    'decreases^farnesylation': Defarnesylation,
    'increases^palmitoylation': Palmitoylation,
    'decreases^palmitoylation': Depalmitoylation,
    'increases^ribosylation': Ribosylation,
    'decreases^ribosylation': Deribosylation,
    'increases^sumoylation': Sumoylation,
    'decreases^sumoylation': Desumoylation,
    # For gene/disease and chemical/disease effects
    'therapeutic': Inhibition
}


class CTDProcessor:
    pass


class CTDChemicalDiseaseProcessor(CTDProcessor):
    """Processes chemical-disease relationships from CTD."""
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        df = self.df[self.df[5] != '']
        for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
            chem_name, chem_mesh_id, chem_cas_id, disease_name, disease_id,\
                direct_ev, inf_gene, inf_score, omim_ids, pmids = list(row)
            if not direct_ev:
                continue
            chem_agent = get_chemical_agent(chem_name, chem_mesh_id,
                                            chem_cas_id)
            disease_agent = get_disease_agent(disease_name, disease_id)
            stmt_types = get_statement_types(direct_ev)
            for rel_str, stmt_type in stmt_types.items():
                anns = {'direct_evidence': rel_str}
                evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns)
                       for pmid in pmids.split('|')]
                stmt = stmt_type(chem_agent, disease_agent,
                                 evidence=evs)
                self.statements.append(stmt)


class CTDGeneDiseaseProcessor(CTDProcessor):
    """Processes gene-disease relationships from CTD."""
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        df = self.df[self.df[4] != '']
        for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
            gene_name, gene_entrez_id, disease_name, disease_id, direct_ev, \
                inf_chem, inf_score, omim_ids, pmids = list(row)
            if not direct_ev:
                continue
            disease_agent = get_disease_agent(disease_name, disease_id)
            gene_agent = get_gene_agent(gene_name, gene_entrez_id)
            stmt_types = get_statement_types(direct_ev)
            for rel_str, stmt_type in stmt_types.items():
                anns = {'direct_evidence': rel_str}
                evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns)
                       for pmid in pmids.split('|')]
                stmt = stmt_type(gene_agent, disease_agent,
                                 evidence=evs)
                self.statements.append(stmt)


class CTDChemicalGeneProcessor(CTDProcessor):
    """Processes chemical-gene relationships from CTD."""
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        for _, row in tqdm.tqdm(self.df.iterrows(), total=len(self.df)):
            chem_name, chem_mesh_id, chem_cas_id, gene_name, gene_entrez_id, \
                gene_forms, organism_name, organism_tax_id, txt, \
                rels, pmids = list(row)

            chem_agent = get_chemical_agent(chem_name, chem_mesh_id,
                                            chem_cas_id)
            gene_agent = get_gene_agent(gene_name, gene_entrez_id)
            stmt_types = get_statement_types(rels)
            context = get_context(organism_name, organism_tax_id)
            for rel_str, stmt_type in stmt_types.items():
                anns = {'interaction_action': rel_str}
                evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns,
                                context=context)
                       for pmid in pmids.split('|')]
                stmt = stmt_type(chem_agent, gene_agent, evidence=evs)
                self.statements.append(stmt)


def get_context(organism_name, organism_tax_id):
    if not organism_tax_id:
        return None
    tax_id = str(int(organism_tax_id))
    return RefContext(organism_name,
                      db_refs={'TAXONOMY': tax_id})


def get_statement_types(rel_str):
    rels = rel_str.split('|')
    return {rel: rel_mapping[rel] for rel in rels if rel in rel_mapping}


def get_disease_agent(name, disease_id):
    groundings = disease_id.split('|')
    db_refs = {}
    for gr in groundings:
        db_ns, db_id = gr.split(':')
        db_refs[db_ns] = db_id
    standard_name, db_refs = standardize_name_db_refs(db_refs)
    if standard_name:
        name = standard_name
    return Agent(name, db_refs=db_refs)


def get_gene_agent(name, gene_entrez_id):
    db_refs = {'EGID': gene_entrez_id}
    hgnc_id = hgnc_client.get_hgnc_id(name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
    standard_name, db_refs = standardize_name_db_refs(db_refs)
    if standard_name:
        name = standard_name
    return Agent(name, db_refs=db_refs)


def get_chemical_agent(name, mesh_id, cas_id):
    db_refs = {'MESH': mesh_id}
    if cas_id:
        db_refs['CAS'] = cas_id
    db_refs = standardize_db_refs(db_refs)
    return Agent(name, db_refs=db_refs)