import tqdm
from indra.statements import *
from indra.databases import hgnc_client
from indra.statements.validate import assert_valid_db_refs
from indra.ontology.standardize import standardize_db_refs, get_standard_agent


# These mappings are only relevant for chemical-gene relations, for
# gene-disease and chemical-disease relations, the only two types are
# therapeutic (mapped to Inhibition) and marker/mechanism (unmapped).
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
}


class CTDProcessor:
    """Parent class for CTD relation-specific processors."""
    def __init__(self, df):
        self.df = df
        self.statements = []


class CTDChemicalDiseaseProcessor(CTDProcessor):
    """Processes chemical-disease relationships from CTD."""

    def extract_statements(self):
        df = self.df[self.df[5] == 'therapeutic']
        for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
            chem_name, chem_mesh_id, chem_cas_id, disease_name, disease_id,\
                direct_ev, inf_gene, inf_score, omim_ids, pmids = list(row)
            if not direct_ev:
                continue
            chem_agent = get_chemical_agent(chem_name, chem_mesh_id,
                                            chem_cas_id)
            disease_agent = get_disease_agent(disease_name, disease_id)
            anns = {'direct_evidence': 'therapeutic'}
            evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns)
                   for pmid in pmids.split('|')]
            stmt = Inhibition(chem_agent, disease_agent,
                              evidence=evs)
            self.statements.append(stmt)


class CTDGeneDiseaseProcessor(CTDProcessor):
    """Processes gene-disease relationships from CTD."""

    def extract_statements(self):
        df = self.df[self.df[4] == 'therapeutic']
        for _, row in tqdm.tqdm(df.iterrows(), total=len(df)):
            gene_name, gene_entrez_id, disease_name, disease_id, direct_ev, \
                inf_chem, inf_score, omim_ids, pmids = list(row)
            if not direct_ev:
                continue
            disease_agent = get_disease_agent(disease_name, disease_id)
            gene_agent = get_gene_agent(gene_name, gene_entrez_id)
            anns = {'direct_evidence': 'therapeutic'}
            evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns)
                   for pmid in pmids.split('|')]
            stmt = Inhibition(gene_agent, disease_agent,
                              evidence=evs)
            self.statements.append(stmt)


class CTDChemicalGeneProcessor(CTDProcessor):
    """Processes chemical-gene relationships from CTD."""

    def extract_statements(self):
        for _, row in tqdm.tqdm(self.df.iterrows(), total=len(self.df)):
            chem_name, chem_mesh_id, chem_cas_id, gene_name, gene_entrez_id, \
                gene_forms, organism_name, organism_tax_id, txt, \
                rels, pmids = list(row)

            chem_agent = get_chemical_agent(chem_name, chem_mesh_id,
                                            chem_cas_id)
            gene_agent = get_gene_agent(gene_name, gene_entrez_id)
            stmt_types = self.get_statement_types(rels, chem_name, txt)
            context = get_context(organism_name, organism_tax_id)
            for rel_str, stmt_type in stmt_types.items():
                anns = {'interaction_action': rel_str}
                evs = [Evidence(source_api='ctd', pmid=pmid, annotations=anns,
                                context=context)
                       for pmid in pmids.split('|')]
                stmt = stmt_type(chem_agent, gene_agent, evidence=evs)
                self.statements.append(stmt)

    @staticmethod
    def get_statement_types(rel_str, chem_name, txt):
        rels = rel_str.split('|')
        reactions = {rel for rel in rels if 'reaction' in rel}
        # If there is a reaction involved and the chemical is not the first
        # element of the reaction description then it is embedded and should
        # not be picked up here (since when there are patterns like
        # A->[B->C], we can pick up B->C from its own separate row).
        if reactions and not txt.startswith(chem_name):
            return {}

        # We now map the relations to INDRA Statements
        mapped_rels = {rel: rel_mapping[rel]
                       for rel in rels if rel in rel_mapping}

        # If we have a decreases^reaction and we know that the chemical name
        # is the first one in the description then we have something like
        # A-|[B->C] meaning that we will have to flip the polarity to
        # capture the fact that the reaction is decreased.
        if 'decreases^reaction' in reactions:
            mapped_rels = {rel: get_inverse_stmt(stmt_type)
                           for rel, stmt_type in mapped_rels.items()}
        return mapped_rels


def get_inverse_stmt(stmt_type):
    if issubclass(stmt_type, Modification):
        return modclass_to_inverse[stmt_type]
    elif stmt_type == Activation:
        return Inhibition
    elif stmt_type == Inhibition:
        return Activation
    elif stmt_type == IncreaseAmount:
        return DecreaseAmount
    elif stmt_type == DecreaseAmount:
        return IncreaseAmount
    raise ValueError('Unexpected statement type %s' % stmt_type)


def get_context(organism_name, organism_tax_id):
    if not organism_tax_id:
        return None
    tax_id = str(int(organism_tax_id))
    db_refs = {'TAXONOMY': tax_id}
    assert_valid_db_refs(db_refs)
    species = RefContext(organism_name,
                         db_refs=db_refs)
    bc = BioContext(species=species)
    return bc


def get_disease_agent(name, disease_id):
    groundings = disease_id.split('|')
    db_refs = {}
    for gr in groundings:
        db_ns, db_id = gr.split(':')
        db_refs[db_ns] = db_id
    return get_standard_agent(name, db_refs)


def get_gene_agent(name, gene_entrez_id):
    db_refs = {'EGID': gene_entrez_id}
    hgnc_id = hgnc_client.get_hgnc_id(name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
    return get_standard_agent(name, db_refs)


def get_chemical_agent(name, mesh_id, cas_id):
    db_refs = {'MESH': mesh_id}
    if cas_id:
        db_refs['CAS'] = cas_id
    return get_standard_agent(name, db_refs)