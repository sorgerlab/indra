from indra.statements import *
from indra.databases.hgnc_client import get_hgnc_id
from indra.ontology.standardize import get_standard_agent

rel_mapping = {
    'CATALYSIS': Activation,
    'INHIBITION': Inhibition,
    'HETERODIMER_ASSOCIATION': Complex,
    'CATALYSIS;HETERODIMER_ASSOCIATION': Complex
}


class AcsnProcessor:
    def __init__(self, relations_df, correspondence_dict,
                 fplx_lookup):
        self.relations_df = relations_df
        self.correspondence_dict = correspondence_dict
        self.fplx_lookup = fplx_lookup
        self.statements = []

    def extract_statements(self):
        # This is where we implement the Statement extraction
        for _, row in self.relations_df.iterrows():
            agA, stmt_types, agB, pmids = list(row)

            agA = get_agent(agA, self.correspondence_dict,
                            self.fplx_lookup)
            agB = get_agent(agB, self.correspondence_dict,
                            self.fplx_lookup)
            stmt_type = get_stmt_type(stmt_types)
            if stmt_type:
                if agA and agB:
                    if str(pmids) == 'nan':
                        evs = [Evidence(source_api='acsn')]

                    else:
                        evs = [Evidence(source_api='acsn', pmid=pmid)
                               for pmid in pmids.split(';')]

                    if stmt_type == Complex:
                        stmt = stmt_type([agA, agB], evidence=evs)
                    else:
                        stmt = stmt_type(agA, agB, evidence=evs)

                    self.statements.append(stmt)


def get_agent(ag_name, correspondence_dict, fplx_lookup):
    db_refs = {}
    if ag_name in correspondence_dict and len(correspondence_dict[ag_name]) == 1:
        grounded_gene = next(iter(correspondence_dict[ag_name]))
        hgnc_id = get_hgnc_id(grounded_gene)
        if hgnc_id:
            db_refs['HGNC'] = hgnc_id
        db_refs['TEXT'] = grounded_gene
        return Agent(grounded_gene, db_refs=db_refs)

    elif ag_name in correspondence_dict and len(correspondence_dict) > 1:
        fplx_rel = fplx_lookup.get(tuple(sorted(correspondence_dict[ag_name])))
        if fplx_rel:
            db_refs['FPLX'] = fplx_rel
            return Agent(fplx_rel, db_refs=db_refs)


def get_stmt_type(stmt_type):
    for r in rel_mapping:
        if stmt_type in rel_mapping:
            return rel_mapping[stmt_type]
