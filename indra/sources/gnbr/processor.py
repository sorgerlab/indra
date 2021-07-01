"""This module contains the processor for GNBR. There are several, each
corresponding to different kinds of interactions."""
import re
import pandas as pd
from indra.statements import *
from indra.ontology.standardize import get_standard_agent


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


cheby_pattern = re.compile(r'^CHEBI:(\d+)$')
cheby_no_prefix_pattern = re.compile(r'^\d+$')

mesh_pattern = re.compile(r'^MESH:([CD]\d+)$')
mesh_no_prefix_pattern = re.compile(r'^[CD]\d+$')

entrez_pattern = re.compile(r'^(\d+)$')
entrez_with_tax_pattern = re.compile(r'^(\d+)\(Tax:(\d+)\)$')

omim_pattern = re.compile(r'^OMIM:(\d+)$')
omim_no_prefix_pattern = re.compile(r'^(\d+)$')


class GnbrProcessor:
    """A processor for interactions in the GNBR dataset.

    Parameters
    ----------
    df1 :
        Dataframe of dependency paths and themes.
    df2 :
        Dataframe of dependency paths and agents.
    first_type :
        The type of the first entity in the data frame.
    second_type :
        The type of the second entity in the data frame.
    """
    def __init__(self, df1: pd.DataFrame, df2: pd.DataFrame,
                 first_type: str, second_type: str,
                 indicator_only: bool = True) -> None:
        self.df1 = df1
        self.df2 = df2
        self.df2.columns = ['id', 'sentence_num', 'nm_1_form', 'nm_1_loc',
                            'nm_2_form', 'nm_2_loc', 'nm_1_raw', 'nm_2_raw',
                            'nm_1_dbid', 'nm_2_dbid', '1_type', '2_type',
                            'path', 'sentence']
        self.df2['path'] = df2['path'].str.lower()
        self.first_type = first_type
        self.second_type = second_type
        self.indicator_only = indicator_only
        self.statements = []

    def extract_stmts(self):
        """Extend the statements list with mappings."""
        if self.first_type == 'gene' and self.second_type == 'gene':
            statement_mappings = gene_gene_stmt_mappings
        elif self.first_type == 'chemical' and self.second_type == 'gene':
            statement_mappings = chem_gene_stmt_mappings
        elif self.first_type == 'gene' and self.second_type == 'disease':
            statement_mappings = gene_disease_stmt_mappings
        else:
            statement_mappings = chem_disease_stmt_mappings
        for rel_type, stmt_type in statement_mappings.items():
            constraint = (self.df1[rel_type] > 0)
            if self.indicator_only:
                constraint &= (self.df1['%s.ind' % rel_type] == 1)
            df_part = self.df1[constraint]
            self.statements.extend(self._extract_stmts_by_class(df_part,
                                                                stmt_type))

    def _extract_stmts_by_class(self, df, stmt_class):
        """Make a given class of Statements from a subset of the dataframe.

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
            else:
                agent1 = get_std_chemical(row['nm_1_raw'], row['nm_1_dbid'])

            if self.second_type == 'gene':
                agent2 = get_std_gene(row['nm_2_raw'], row['nm_2_dbid'])
            else:
                agent2 = get_std_disease(row['nm_2_raw'], row['nm_2_dbid'])

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
    :
        A standardized Agent object.
    """
    db_refs = {'TEXT': raw_string} if not pd.isna(raw_string) else {}
    name = raw_string if not pd.isna(raw_string) else db_id
    # TODO: Split into multiple agents if ; is encountered
    if not pd.isna(db_id) and ';' in db_id:
        db_id = db_id.split(';')[0]

    if pd.isna(db_id):
        pass
    elif entrez_pattern.match(db_id):
        db_refs['EGID'] = db_id
    else:
        match = entrez_with_tax_pattern.match(db_id)
        if not match:
            raise ValueError('Unexpected gene identifier: %s' % db_id)
        db_refs['EGID'] = match.groups()[0]
    return get_standard_agent(name, db_refs)


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
    :
        A standardized Agent object.
    """
    db_refs = {'TEXT': raw_string} if not pd.isna(raw_string) else {}
    name = raw_string if not pd.isna(raw_string) else db_id
    # TODO: the interpretation of multiple IDs separated by | should
    # be clarified
    if not pd.isna(db_id) and '|' in db_id:
        db_id = db_id.split('|')[0]

    if pd.isna(db_id):
        pass
    elif cheby_pattern.match(db_id):
        db_refs['CHEBI'] = db_id[6:]
    elif cheby_no_prefix_pattern.match(db_id):
        db_refs['CHEBI'] = db_id
    elif mesh_pattern.match(db_id):
        db_refs['MESH'] = db_id[5:]
    elif mesh_no_prefix_pattern.match(db_id):
        db_refs['MESH'] = db_id
    else:
        raise ValueError('Unexpected chemical identifier: %s' % db_id)
    return get_standard_agent(name, db_refs)


def get_std_disease(raw_string: str, db_id: str) -> Agent:
    """Standardize disease names.

    Parameters
    ----------
    raw_string :
        Name of the agent in the GNBR dataset.
    db_id :
        Entrez identifier of the agent.

    Returns
    -------
    :
        A standardized Agent object.
    """
    db_refs = {'TEXT': raw_string} if not pd.isna(raw_string) else {}
    name = raw_string if not pd.isna(raw_string) else db_id

    if pd.isna(db_id):
        pass
    elif omim_no_prefix_pattern.match(db_id):
        db_refs['OMIM'] = db_id
    elif omim_pattern.match(db_id):
        db_refs['OMIM'] = db_id[5:]
    elif mesh_no_prefix_pattern.match(db_id):
        db_refs['MESH'] = db_id
    elif mesh_pattern.match(db_id):
        db_refs['MESH'] = db_id[5:]
    else:
        raise ValueError('Unexpected disease identifier: %s' % db_id)
    return get_standard_agent(name, db_refs)


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
