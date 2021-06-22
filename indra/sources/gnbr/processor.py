from indra.statements import *
from indra.statements import Agent
from indra.ontology.standardize import standardize_agent_name


class GnbrGeneGeneProcessor:
    def __init__(self, df1, df2):
        self.df1 = df1
        self.df2 = df2
        self.df2.columns = ['id', 'sentence_num', 'nm_1_form', 'nm_1_loc',
                            'nm_2_form', 'nm_2_loc', 'nm_1_raw', 'nm_2_raw',
                            'nm_1_dbid', 'nm_2_dbid', '1_type', '2_type',
                            'path', 'sentence']
        self.df2['path'] = df2['path'].str.lower()
        self.statements = []

    def extract_activations(self):
        df1_activations = self.df1[(self.df1['V+.ind'] == 1) & (self.df1['V+'] > 0)]
        df = df1_activations.join(self.df2.set_index('path'), on='path')

        for index, row in df.iterrows():
            agent1 = self.standardize_agent(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.standardize_agent(row['nm_2_raw'], row['nm_2_dbid'])
            self.statements.append(Activation(agent1, agent2))

    @staticmethod
    def standardize_agent(raw_string, db_id):
        agent = Agent(raw_string, db_refs={'EGID': db_id, 'TEXT': raw_string})
        standardize_agent_name(agent)
        return agent