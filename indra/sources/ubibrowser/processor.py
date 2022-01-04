from indra.statements import *
from indra.ontology.standardize import get_standard_agent


class UbiBrowserProcessor:
    """Processor for UbiBrowser data."""
    def __init__(self, e3_df, dub_df):
        self.e3_df = e3_df
        self.dub_df = dub_df
        self.statements = []

    def extract_statements(self):
        for df, stmt_type in [(self.e3_df, Ubiquitination),
                              (self.dub_df, Deubiquitination)]:
            for _, row in df.iterrows():
                stmt = self._process_row(row, stmt_type)
                if stmt:
                    self.statements.append(stmt)

    @staticmethod
    def _process_row(row, stmt_type):
        # Note that even in the DUB table the subject of the statement
        # is called "E3"
        # There are some examples where a complex is implied (e.g., BMI1-RNF2),
        # for simplicity we just ignore these
        if '-' in row['E3AC']:
            return None
        subj_agent = get_standard_agent(row['E3GENE'], {'UP': row['E3AC']})
        obj_agent = get_standard_agent(row['SUBGENE'], {'UP': row['SUBAC']})
        if row['SOURCE'] == 'MEDLINE' and row['SOURCEID'] != 'UNIPROT':
            # Note: we sometimes get int here
            pmid = str(row['SOURCEID'])
            text = row['SENTENCE']
        else:
            pmid = None
            text = None
        ev = Evidence(source_api='ubibrowser', pmid=pmid, text=text)
        stmt = stmt_type(subj_agent, obj_agent, evidence=[ev])
        return stmt
