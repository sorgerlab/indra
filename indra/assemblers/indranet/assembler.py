import logging
import pandas as pd
from .net import IndraNet
from indra.statements import *
from itertools import permutations
from collections import OrderedDict


logger = logging.getLogger(__name__)
NS_PRIORITY_LIST = (
    'FPLX', 'HGNC', 'UP', 'CHEBI', 'GO', 'MESH', 'HMDB', 'PUBCHEM')


class IndraNetAssembler():
    """Assembler to create an IndraNet object from a list of INDRA statements.

    Parameters
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.

    Attributes
    ----------
    model : IndraNet
        An IndraNet graph object assembled by this class.
    """
    def __init__(self, statements=None):
        self.statements = statements if statements else []
        self.model = None

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        for stmt in stmts:
            self.statements.append(stmt)

    def make_model(self, exclude_stmts=None, complex_members=3):
        """Assemble an IndraNet graph object.

        Parameters
        ----------
        exclude_stmts : list[str]
            A list of statement type names to not include in the graph.
        complex_members : int
            Maximum allowed size of a complex to be included in the graph.
            All complexes larger than complex_members will be rejected. For
            accepted complexes, all permutations of their members will be added
            as edges.

        Returns
        -------
        model : IndraNet
            IndraNet graph object.
        """
        df = self.make_df(exclude_stmts, complex_members)
        model = IndraNet.from_df(df)
        return model

    def make_df(self, exclude_stmts=None, complex_members=3):
        """Create a dataframe containing information extracted from assembler's
        list of statements necessary to build an IndraNet.

        Parameters
        ----------
        exclude_stmts : list[str]
            A list of statement type names to not include into a dataframe.
        complex_members : int
            Maximum allowed size of a complex to be included in the data
            frame. All complexes larger than complex_members will be rejected.
            For accepted complexes, all permutations of their members will be
            added as dataframe records.

        Returns
        -------
        df : pd.DataFrame
            Pandas DataFrame object containing information extracted from
            statements.
        """
        rows = []
        if exclude_stmts:
            exclude_types = tuple(
                get_statement_by_name(st_type) for st_type in exclude_stmts)
        else:
            exclude_types = ()
        for stmt in self.statements:
            # Exclude statements from given exclude list
            if isinstance(stmt, exclude_types):
                logger.warning('Skipping a statement of a type %s.'
                               % type(stmt).__name__)
                continue
            agents = stmt.agent_list()
            # Exclude statements with less than 2 agents
            if len(agents) < 2:
                logger.warning('Skipping a statement with less than 2 agents.')
                continue
            # Handle complexes
            if isinstance(stmt, Complex):
                # Do not add complexes with more members than complex_members
                if len(agents) > complex_members:
                    logger.warning('Skipping a complex with %d members.'
                                   % len(agents))
                    continue
                else:
                    # add every permutation
                    pairs = permutations(agents, 2)
            else:
                pairs = [agents]
            for (agA, agB) in pairs:
                if agA is None or agB is None:
                    logger.warning(
                        'Skipping a statement because agent is None.')
                    continue

                def get_ag_ns_id(ag):
                    # get one pair of ns-id from agent db_refs
                    ag_ns, ag_id = None, None
                    for ns, id in ag.db_refs.items():
                        if ns in NS_PRIORITY_LIST:
                            ag_ns, ag_id = ns, id
                            break
                        if not ns:
                            ag_ns = 'TEXT'
                            ag_id = agent.name
                    return ag_ns, ag_id

                agA_ns, agA_id = get_ag_ns_id(agA)
                agB_ns, agB_id = get_ag_ns_id(agB)
                stmt_type = type(stmt).__name__
                row = OrderedDict([
                    ('agA_name', agA.name),
                    ('agB_name', agB.name),
                    ('agA_ns', agA_ns),
                    ('agA_id', agA_id),
                    ('agB_ns', agB_ns),
                    ('agB_id', agB_id),
                    ('stmt_type', stmt_type),
                    ('evidence_count', len(stmt.evidence)),
                    ('hash', stmt.get_hash(refresh=True)),
                    ('belief', stmt.belief)])
                rows.append(row)
        df = pd.DataFrame.from_dict(rows)
        return df
