from indra.statements import *
from indra.ontology.standardize import get_standard_agent


class UbiBrowserProcessor:
    """Processor for UbiBrowser data."""
    def __init__(self, e3_df, dub_df):
        self.e3_df = e3_df
        self.dub_df = dub_df
        self.statements = []

    def extract_statements(self):
        for df, stmt_type, subj_suffix in \
                [(self.e3_df, Ubiquitination, 'E3'),
                 (self.dub_df, Deubiquitination, 'DUB')]:
            for _, row in df.iterrows():
                stmt = self._process_row(row, stmt_type, subj_suffix)
                if stmt:
                    self.statements.append(stmt)

    @staticmethod
    def _process_row(row, stmt_type, subj_suffix):
        # Note that even in the DUB table the subject of the statement
        # is called "E3"
        # There are some examples where a complex is implied (e.g., BMI1-RNF2),
        # for simplicity we just ignore these
        if '#' in row[f'SwissProt AC ({subj_suffix})']:
            return None
        # Interestingly, some of the E3s are missing entirely, we skip these
        elif row[f'SwissProt AC ({subj_suffix})'] == '-':
            return None
        # Some of the same corner cases apply to the substrate as well
        if row['SwissProt AC (Substrate)'] == '-':
            return None
        subj_agent = \
            get_standard_agent(row[f'Gene Symbol ({subj_suffix})'],
                               {'UP': row[f'SwissProt AC ({subj_suffix})']})
        obj_agent = get_standard_agent(row['Gene Symbol (Substrate)'],
                                       {'UP': row['SwissProt AC (Substrate)']})
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
