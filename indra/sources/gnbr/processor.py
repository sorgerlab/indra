"""This module contains the processor for GNBR. There are several, each
corresponding to different kinds of interactions."""

from indra.statements import *
from indra.statements import Agent
from indra.statements import Evidence
from indra.ontology.standardize import standardize_agent_name
import pandas as pd
import re


gene_gene_stmt_mappings = {
    'V+': Activation,
    'E+': IncreaseAmount,
    'Q':  IncreaseAmount,
    'H':  Complex
}

chem_gene_stmt_mappings = {
    'A+': Activation,
    'A-': Inhibition,
    'N':  Inhibition,
    'B':  Complex,
    'E-': DecreaseAmount
}

gene_disease_stmt_mappings = {
    'Te': Inhibition,
    'G':  Activation
}

chem_disease_stmt_mappings = {
    'T':  Inhibition,
    'C':  Inhibition,
    'Pr': Inhibition,
    'Pa': Inhibition
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
        """Extend the statements list with mappings."""
        if self.first_type == 'gene' and self.second_type == 'gene':
            statement_mappings = gene_gene_stmt_mappings
        elif self.first_type == 'chemical' and self.second_type == 'gene':
            statement_mappings = chem_gene_stmt_mappings
        for rel_type, stmt_type in statement_mappings.items():
            df_part = self.df1[(self.df1['%s.ind' % rel_type] == 1) &
                               (self.df1[rel_type] > 0)]
            self.statements.extend(self._extract_stmts(df_part, stmt_type))

    def _extract_stmts(self, df, stmt_class):
        """Make Statements from the dataframes.

        Parameters
        ----------
        df :
            Filtered dataframe to one particular relationship theme.
        stmt_class :
            Statement type matched to the type of the filtered dataframe.

        Yields
        ------
        stmt :
            Statements produced from the dataframes.
        """
        df_joint = df.join(self.df2.set_index('path'), on='path')
        for index, row in df_joint.iterrows():
            agent1: Agent
            agent2: Agent
            if self.first_type == 'gene':
                agent1 = get_std_gene(row['nm_1_raw'], row['nm_1_dbid'])
            elif self.first_type == 'chemical':
                agent1 = get_std_chemical(row['nm_1_raw'], row['nm_1_dbid'])
            # elif self.first_type == 'disease':
                # agent1 = get_std_disease(row['nm_1_raw'], row['nm_1_dbid'])

            if self.second_type == 'gene':
                agent2 = get_std_gene(row['nm_2_raw'], row['nm_2_dbid'])
            elif self.second_type == 'chemical':
                agent2 = get_std_chemical(row['nm_2_raw'], row['nm_2_dbid'])
            # elif self.second_type == 'disease':
                # agent2 = get_std_disease(row['nm_2_raw'], row['nm_2_dbid'])

            evidence = get_evidence(row)
            if stmt_class == Complex:
                stmt = Complex([agent1, agent2], evidence=evidence)
            else:
                stmt = stmt_class(agent1, agent2, evidence=evidence)
            yield stmt


def get_std_gene(raw_string: str, db_id: str) -> Agent:
    """Standardize gene names.

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
    if re.match('^(\d+)$', db_id):
        match = re.match('^(\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'EGID': match.groups()[0],
                                           'TEXT': raw_string})
    else:
        match = re.match('^(\d+)\(Tax:(\d+)\)$', db_id)
        agent = Agent(raw_string, db_refs={'EGID': match.groups()[0],
                                           'TEXT': raw_string})
    standardize_agent_name(agent)
    return agent


def get_std_chemical(raw_string: str, db_id: str) -> Agent:
    """Standardize chemical names.

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
    if re.match('^CHEBI:(\d+)$', db_id):
        match = re.match('^CHEBI:(\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'CHEBI': match.groups()[0],
                                           'TEXT': raw_string})
    elif re.match('^MESH:([A-Z]\d+)$', db_id):
        match = re.match('^MESH:([A-Z]\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'MESH': match.groups()[0],
                                           'TEXT': raw_string})
    elif db_id == 'nan':
        agent = Agent(name=raw_string, db_refs={'TEXT': raw_string})

    return agent


def get_std_disease(raw_string: str, db_id: str):
    """Standardize disease names.

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
    if re.match('^(\d+)$', db_id):
        match = re.match('^(\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'OMIM': match.groups()[0],
                                           'TEXT': raw_string})
        standardize_agent_name(agent)
    elif re.match('^OMIM:(\d+)$', db_id):
        match = re.match('^OMIM:(\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'OMIM': match.groups()[0],
                                           'TEXT': raw_string})
    elif re.match('^([A-Z]\d+)$', db_id):
        match = re.match('^([A-Z]\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'MESH': match.groups()[0],
                                           'TEXT': raw_string})
    else:
        match = re.match('^MESH:([A-Z]\d+)$', db_id)
        agent = Agent(raw_string, db_refs={'MESH': match.groups()[0],
                                           'TEXT': raw_string})


def get_evidence(row):
    """Give evidence for a Statement.

    Parameters
    ----------
    row :
        Currently investigated row of the dataframe.

    Returns
    -------
    evidence :
        Evidence object with the source_api, the PMID and the original
        sentence.
    """
    pmid = str(row['id']) if row['id'] else None
    evidence = Evidence(source_api='gnbr',
                        pmid=pmid,
                        text=row['sentence'],
                        text_refs={'PMID': pmid})
    return evidence
