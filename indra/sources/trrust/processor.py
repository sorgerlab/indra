from copy import deepcopy
from indra.databases import hgnc_client
from indra.statements import Agent, IncreaseAmount, DecreaseAmount, Evidence


class TrrustProcessor(object):
    """Processor to extract INDRA Statements from Trrust data frame.

    Attributes
    ----------
    df : pandas.DataFrame
        The Trrust table to process.
    statements : list[indra.statements.Statement]
        The list of INDRA Statements extracted from the table.
    """
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        """Process the table to extract Statements."""
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
    """Return a Statement based on its type, agents, and PMID."""
    ev = Evidence(source_api='trrust', pmid=pmid)
    return stmt_cls(deepcopy(tf_agent), deepcopy(target_agent),
                    evidence=[ev])


def get_grounded_agent(gene_name):
    """Return a grounded Agent based on an HGNC symbol."""
    db_refs = {'TEXT': gene_name}
    if gene_name in hgnc_map:
        gene_name = hgnc_map[gene_name]
    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
    if not hgnc_id:
        hgnc_id = hgnc_client.get_current_hgnc_id(gene_name)
    if hgnc_id:
        db_refs['HGNC'] = hgnc_id
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        if up_id:
            db_refs['UP'] = up_id
    agent = Agent(gene_name, db_refs=db_refs)
    return agent


hgnc_map = {
    'CTGF': 'CCN2',
    'CYR61': 'CCN1',
    'MKL1': 'MRTFA',
    'NOV': 'CCN3',
    'RFWD2': 'COP1',
    'SALL4A': 'SALL4',
    'STAT5': 'STAT5A',
    'TRAP': 'ACP5',
    'AES': 'TLE5',
    'SEPT7': 'SEPTIN7'
}
