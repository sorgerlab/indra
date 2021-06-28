"""This module contains the processor for GNBR. There are several, each
corresponding to different kinds of interactions."""

from indra.statements import *
from indra.statements import Agent
from indra.statements import Evidence
from indra.ontology.standardize import standardize_agent_name
import pandas as pd


statement_mappings = {
    'V+': Activation,
    'E+': IncreaseAmount,
    'Q':  IncreaseAmount,
    'H':  Complex,
    'A+': Activation,
    'A-': Inhibition,
    'N':  Inhibition,
    'B':  Complex,
    'E-': DecreaseAmount,
}


class GnbrProcessor:
    """A processor for interactions in the GNBR dataset.

    Parameters
    ----------
    df1 :
        Dataframe of dependency paths and themes.
    df2 :
        Dataframe of dependency paths and agents.
    """
    def __init__(self, df1: pd.DataFrame, df2: pd.DataFrame,
                 first_type: str, second_type: str) -> None:
        self.df1 = df1
        self.df2 = df2
        self.df2.columns = ['id', 'sentence_num', 'nm_1_form', 'nm_1_loc',
                            'nm_2_form', 'nm_2_loc', 'nm_1_raw', 'nm_2_raw',
                            'nm_1_dbid', 'nm_2_dbid', '1_type', '2_type',
                            'path', 'sentence']
        self.df2['path'] = df2['path'].str.lower()
        self.first_type = first_type
        self.second_type = second_type
        self.statements = []

    def extract_stmts(self):
        for rel_type, stmt_type in statement_mappings.items():
            df_part = self.df1[(self.df1['%s.ind' % rel_type] == 1) &
                               (self.df1[rel_type] > 0)]
            self.statements.extend(self._extract_stmts(df_part, stmt_type))

    def _extract_stmts(self, df, stmt_class):
        df_joint = df.join(self.df2.set_index('path'), on='path')
        for index, row in df_joint.iterrows():
            if self.first_type == 'gene':
                agent1 = get_std_gene(row['nm_1_raw'], row['nm_1_dbid'])
            elif self.first_type == 'chemical':
                agent1 = get_std_chemical(row['nm_1_raw'], row['nm_1_dbid'])
            elif self.first_type == 'disease':
                agent1 = get_std_disease(row['nm_1_raw'], row['nm_1_dbid'])

            if self.second_type == 'gene':
                agent2 = get_std_gene(row['nm_2_raw'], row['nm_2_dbid'])
            elif self.second_type == 'chemical':
                agent2 = get_std_chemical(row['nm_2_raw'], row['nm_2_dbid'])
            elif self.second_type == 'disease':
                agent2 = get_std_disease(row['nm_2_raw'], row['nm_2_dbid'])

            evidence = get_evidence(row)
            if stmt_class == Complex:
                stmt = Complex([agent1, agent2], evidence=evidence)
            else:
                stmt = stmt_class(agent1, agent2, evidence=evidence)
            yield stmt


def get_std_gene(raw_string: str, db_id: str) -> Agent:
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
    agent: Agent = Agent(raw_string, db_refs={'EGID': db_id,
                                              'TEXT': raw_string})
    standardize_agent_name(agent)
    return agent


def get_std_chemical(raw_string: str, db_id: str) -> Agent:
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
        agent = Agent(name=raw_string, db_refs={'TEXT': raw_string})
    elif db_id.startswith('CHEBI:'):
        agent = Agent(raw_string, db_refs={'CHEBI': db_id,
                                           'TEXT': raw_string})
        standardize_agent_name(agent)
    elif db_id.startswith('MESH:'):
        agent = Agent(raw_string, db_refs={'MESH': db_id.split(':')[1],
                                           'TEXT': raw_string})
        standardize_agent_name(agent)
    elif '(Tax:' in db_id:
        agent = Agent(raw_string, db_refs={'MESH': db_id.split('(')[0],
                                           'TEXT': raw_string})
        standardize_agent_name(agent)
    else:
        agent = Agent(raw_string, db_refs={'MESH': db_id,
                                           'TEXT': raw_string})
        standardize_agent_name(agent)

    return agent


def get_std_disease(raw_string: str, db_id: str):
    """Standardize agent (chemical) names."""


def get_evidence(row):
    pmid = str(row['id']) if row['id'] else None
    evidence = Evidence(source_api='gnbr',
                        pmid=pmid,
                        text=row['sentence'],
                        text_refs={'PMID': pmid})
    return evidence
