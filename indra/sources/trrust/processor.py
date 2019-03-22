from indra.statements import Agent, IncreaseAmount, DecreaseAmount
from indra.databases import hgnc_client, uniprot_client


class TrrustProcessor(object):
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        for _, (tf, target, effect, refs) in self.df.iterrows():
            tf_agent = get_grounded_agent(tf)
            target_agent = get_grounded_agent(target_agent)
            if effect == 'Activation':
                stmt_cls = IncreaseAmount
            elif effect == 'Repression':
                stmt_cls = DecreaseAmount
            else:
                continue
            stmt = stmt_cls(tf_agent, target_agent)
            self.statements.append(stmt)


def get_grounded_agent(gene_name):
    db_refs = {'TEXT': gene_name}
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        if up_id:
            db_refs['UP'] = up_id
    agent = Agent(gene_name, db_refs=db_refs)
    return agent
