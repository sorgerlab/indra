import re
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

    # There's a column that is always a - character
    assert row['dash'] == '-', row['dash']

    exp_method_id, exp_method_name = parse_psi_mi(row['exp_method'])
    int_type__id, int_type_name = parse_psi_mi(row['int_type'])

    annotations = {
        'exp_method': {'id': exp_method_id, 'name': exp_method_name},
        'int_type': {'id': int_type__id, 'name': int_type_name},
    }

    ev = Evidence(source_api='virhostnet', annotations=annotations)

    stmt = Complex([host_agent, vir_agent], evidence=[ev])
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


def parse_psi_mi(psi_mi_str):
    # Example: psi-mi:"MI:0018"(two hybrid)
    match = re.match(r'psi-mi:"(.+)"\((.+)\)')
    mi_id, name = match.groups()
    return mi_id, name
