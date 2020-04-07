from indra.statements import Agent


class VirhostnetProcessor():
    def __init__(self, df):
        self.df = df
        self.statements = []

    def extract_statements(self):
        for _, row in self.df.iterrows():
            host_up, vir_up, _, _, _, _, host_tax, vir_tax, psi_mi_pa, \
                psi_mi_vhn, vhn_ids, score = row
            stmt = process_row(host_up, vir_up, host_tax, vir_tax, psi_mi_pa,
                               psi_mi_vhn, vhn_ids, score)
            if stmt:
                self.statements.append(stmt)


def process_row(host_up, vir_up, host_tax, vir_tax, psi_mi_pa, psi_mi_vhn,
                vhn_ids, score):
    return None
