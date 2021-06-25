"""This module contains the processors for GNBR. There are several, each corresponding
to different kinds of interactions."""

from indra.statements import *
from indra.statements import Agent
from indra.statements import Evidence
from indra.ontology.standardize import standardize_agent_name
import pandas as pd


class GnbrGeneGeneProcessor:
    """A processor for gene-gene interactions in the GNBR dataset.

    Parameters
    ----------
    df1 :
        Dataframe of dependency paths and themes.
    df2 :
        Dataframe of dependency paths and agents.
    """
    def __init__(self, df1: pd.DataFrame, df2: pd.DataFrame) -> None:
        self.df1 = df1
        self.df2 = df2
        self.df2.columns = ['id', 'sentence_num', 'nm_1_form', 'nm_1_loc',
                            'nm_2_form', 'nm_2_loc', 'nm_1_raw', 'nm_2_raw',
                            'nm_1_dbid', 'nm_2_dbid', '1_type', '2_type',
                            'path', 'sentence']
        self.df2['path'] = df2['path'].str.lower()
        self.statements = []

    def extract_activations(self) -> None:
        """Make Activation Statements from the dataframes."""
        df1_activations = self.df1[(self.df1['V+.ind'] == 1) &
                                   (self.df1['V+'] > 0)]
        df = df1_activations.join(self.df2.set_index('path'), on='path')
        for index, row in df.iterrows():
            agent1 = self.standardize_agent(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.standardize_agent(row['nm_2_raw'], row['nm_2_dbid'])
            evidence = Evidence(source_api='gnbr',
                                pmid=row['id'],
                                text=row['sentence'],
                                text_refs={'PMID': row['id']}
                                )
            self.statements.append(Activation(agent1, agent2, evidence))

    def extract_increase_amount(self) -> None:
        """Make IncreaseAmount Statements from the dataframes."""
        df1_increase_amounts = self.df1[((self.df1['E+.ind'] == 1) &
                                         (self.df1['E+'] > 0)) |
                                        ((self.df1['Q.ind'] == 1) &
                                         (self.df1['Q'] > 0))]
        df = df1_increase_amounts.join(self.df2.set_index('path'), on='path')
        for index, row in df.iterrows():
            agent1 = self.standardize_agent(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.standardize_agent(row['nm_2_raw'], row['nm_2_dbid'])
            evidence = Evidence(source_api='gnbr',
                                pmid=row['id'],
                                text=row['sentence'],
                                text_refs={'PMID': row['id']}
                                )
            self.statements.append(IncreaseAmount(agent1, agent2, evidence))

    def extract_complexes(self):
        """Make Complex Statements from the dataframes."""
        df1_complexes = self.df1[(self.df1['H.ind'] == 1) &
                                        (self.df1['H'] > 0)]
        df = df1_complexes.join(self.df2.set_index('path'), on='path')
        for index, row in df.iterrows():
            agent1 = self.standardize_agent(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.standardize_agent(row['nm_2_raw'], row['nm_2_dbid'])
            evidence = Evidence(source_api='gnbr',
                                pmid=row['id'],
                                text=row['sentence'],
                                text_refs={'PMID': row['id']}
                                )
            self.statements.append(Complex([agent1, agent2], evidence))


    @staticmethod
    def standardize_agent(raw_string: str, db_id: str) -> Agent:
        """Standardize agent names.

        Parameters
        ----------
        raw_string :
            Name of the agent in the GNBR dataset.
        db_id :
            Entrez identifier of the agent.

        Returns
        -------
        agent :
            A standardized Agent object.
        """
        agent = Agent(raw_string, db_refs={'EGID': db_id, 'TEXT': raw_string})
        standardize_agent_name(agent)
        return agent

class GnbrChemicalGeneProcessor:
    """A processor for gene-gene interactions in the GNBR dataset.

    Parameters
    ----------
    df1 :
        Dataframe of dependency paths and themes.
    df2 :
        Dataframe of dependency paths and agents.
    """
    def __init__(self, df1: pd.DataFrame, df2: pd.DataFrame) -> None:
        self.df1 = df1
        self.df2 = df2
        self.df2.columns = ['id', 'sentence_num', 'nm_1_form', 'nm_1_loc',
                            'nm_2_form', 'nm_2_loc', 'nm_1_raw', 'nm_2_raw',
                            'nm_1_dbid', 'nm_2_dbid', '1_type', '2_type',
                            'path', 'sentence']
        self.df2['path'] = df2['path'].str.lower()
        self.statements = []

    def extract_activations(self):
        """Make Activation Statements from the dataframes."""
        df1_activations = self.df1[(self.df1['A+.ind'] == 1) &
                                   (self.df1['A+'] > 0)]
        df = df1_activations.join(self.df2.set_index('path'), on='path')
        for index, row in df.iterrows():
            agent1 = self.std_gene(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.std_gene(row['nm_2_raw'], row['nm_2_dbid'])
            evidence = Evidence(source_api='gnbr',
                                pmid=row['id'],
                                text=row['sentence'],
                                text_refs={'PMID': row['id']}
                                )
            self.statements.append(Activation(agent1, agent2, evidence))

    def extract_inhibition(self):
        """Make Inhibition Statements from the dataframes."""
        df1_inhibitions = self.df1[((self.df1['A-.ind'] == 1) &
                                    (self.df1['A-'] > 0)) |
                                   ((self.df1['N.ind'] == 1) &
                                    (self.df1['N'] > 0))]
        df = df1_inhibitions.join(self.df2.set_index('path'), on='path')
        for index, row in df.iterrows():
            agent1 = self.std_chemical(row['nm_1_raw'], row['nm_1_dbid'])
            agent2 = self.std_gene(row['nm_2_raw'], row['nm_2_dbid'])
            evidence = Evidence(source_api='gnbr',
                                pmid=row['id'],
                                text=row['sentence'],
                                text_refs={'PMID': row['id']}
                                )
            self.statements.append(Inhibition(agent1, agent2, evidence))

    @staticmethod
    def std_gene(raw_string: str, db_id: str) -> Agent:
        """Standardize agent (gene) names.

        Parameters
        ----------
        raw_string :
            Name of the agent in the GNBR dataset.
        db_id :
            Entrez identifier of the agent.

        Returns
        -------
        agent :
            A standardized Agent object.
        """
        agent: Agent = Agent(raw_string, db_refs={'EGID': db_id, 'TEXT': raw_string})
        standardize_agent_name(agent)
        return agent

    @staticmethod
    def std_chemical(raw_string: str, db_id: str) -> Agent:
        """Standardize agent (chemical) names.

        Parameters
        ----------
        raw_string :
            Name of the agent in the GNBR dataset.
        db_id :
            Entrez identifier of the agent.

        Returns
        -------
        agent :
            A standardized Agent object.
        """
        agent: Agent
        if db_id == 'null':
            agent.name = raw_string
            agent.db_refs = raw_string
        elif db_id.startswith('CHEBI:'):
            agent = Agent(raw_string, db_refs={'CHEBI': db_id,
                                               'TEXT': raw_string})
            standardize_agent_name(agent)
        elif db_id.startswith('MESH:'):
            agent = Agent(raw_string, db_refs={'MESH': db_id.split(':')[1],
                                               'TEXT': raw_string})
            standardize_agent_name(agent)
        else:
            agent = Agent(raw_string, db_refs={'MESH': db_id,
                                               'TEXT': raw_string})
            standardize_agent_name(agent)

        return agent
