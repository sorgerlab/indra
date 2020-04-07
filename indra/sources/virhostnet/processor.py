from indra.databases import uniprot_client
from indra.statements import Agent, Complex, Evidence
from indra.preassembler.grounding_mapper import standardize_agent_name


class VirhostnetProcessor():
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        for _, row in self.df.iterrows():
            stmt = process_row(row)
            if stmt:
                self.statements.append(stmt)


def process_row(row):
    host_agent = get_agent_from_grounding(row['host_grounding'])
    vir_agent = get_agent_from_grounding(row['vir_grounding'])
    stmt = Complex([host_agent, vir_agent])
    return stmt


def get_agent_from_grounding(grounding):
    db_ns, db_id = grounding.split(':')
    # Assume UniProt or RefSeq IDs
    assert db_ns in {'uniprotkb', 'refseq', 'ddbj/embl/genbank'}, db_ns
    if db_ns == 'uniprotkb':
        if '-' in db_id:
            up_id, feat_id = db_id.split('-')
            # Assume it's a feature ID
            assert feat_id.startswith('PRO'), feat_id
            db_refs = {'UP': up_id, 'UPPRO': feat_id}
        else:
            db_refs = {'UP': db_id}
    elif db_ns == 'refseq':
        db_refs = {'REFSEQ_PROT': db_id}
    else:
        db_refs = {'GENBANK': db_id}
    agent = Agent(db_id, db_refs=db_refs)
    standardize_agent_name(agent, standardize_refs=True)
    return agent
