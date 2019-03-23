from copy import deepcopy
from indra.databases import hgnc_client, uniprot_client
from indra.statements import Agent, IncreaseAmount, DecreaseAmount, Evidence


class TrrustProcessor(object):
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        for _, (tf, target, effect, refs) in self.df.iterrows():
            tf_agent = get_grounded_agent(tf)
            target_agent = get_grounded_agent(target)
            if effect == 'Activation':
                stmt_cls = IncreaseAmount
            elif effect == 'Repression':
                stmt_cls = DecreaseAmount
            else:
                continue
            pmids = refs.split(';')
            for pmid in pmids:
                stmt = make_stmt(stmt_cls, tf_agent, target_agent, pmid)
                self.statements.append(stmt)


def make_stmt(stmt_cls, tf_agent, target_agent, pmid):
    ev = Evidence(source_api='trrust', pmid=pmid)
    return stmt_cls(deepcopy(tf_agent), deepcopy(target_agent),
                    evidence=[ev])


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
