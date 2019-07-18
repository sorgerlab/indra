import logging
import pandas as pd
from .indra_net import IndraNet
from indra.statements import *
from itertools import permutations
from collections import OrderedDict


logger = logging.getLogger()
NS_PRIORITY_LIST = ('FPLX', 'HGNC', 'GO', 'MESH', 'HMDB', 'CHEBI', 'PUBCHEM')
SIGN_DICT = {'Activation': 0, 'Inhibition': 1, 'IncreaseAmount': 0,
             'DecreaseAmount': 1}


class IndranetAssembler():
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

    def make_model(self, signed=False, exclude_stmts=None, complex_members=3):
        """Assemble an IndraNet graph object.

        Parameters
        ----------
        signed : bool
            Whether the edges of a returned graph should be signed.
        exclude_stmts : list[str]
            A list of statement type names to not include into a graph.
        complex_members : int
            A maximum allowed size of a complex to be included in the graph.
            All complexes larger than complex_members will be rejected. For
            accepted complexes, all permutation of their members will be added
            as edges.

        Returns
        -------
        model : IndraNet
            IndraNet graph object.
        """
        df = self.make_df(signed, exclude_stmts, complex_members)
        model = IndraNet.from_df(df)
        return model

    def make_df(self, signed=False, exclude_stmts=None, complex_members=3):
        """Create a data frame containing information extracted from assembler's
        list of statements necessary to build an IndraNet.

        Parameters
        ----------
        signed : bool
            Whether the data frame should contain 'sign' column.
        exclude_stmts : list[str]
            A list of statement type names to not include into a data frame.
        complex_members : int
            A maximum allowed size of a complex to be included in the data
            frame. All complexes larger than complex_members will be rejected.
            For accepted complexes, all permutation of their members will be
            added as data frame records.

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
                continue
            agents = stmt.agent_list()
            # Exclude statements with less than 2 agents
            if len(agents) < 2:
                continue
            # Handle complexes
            if isinstance(stmt, Complex):
                # Do not add complexes with more members than complex_members
                if len(agents) > complex_members:
                    continue
                else:
                    # add every permutation
                    pairs = permutations(agents, 2)
            else:
                pairs = [agents]
            for (agA, agB) in pairs:
                if agA is None or agB is None:
                    continue

                def get_ag_ns_id(ag):
                    # get one pair of ns-id from agent db_refs
                    ns = None
                    for ns, id in ag.db_refs.items():
                        if ns in NS_PRIORITY_LIST:
                            ns, id = ns, id
                            break
                        if not ns:
                            ns = 'TEXT'
                            id = agent.name
                    return ns, id

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
                if signed:
                    try:
                        sign = SIGN_DICT[stmt_type]
                    except KeyError:
                        logger.warning('Could not find a sign for %s. '
                                       'Using sign=0 by default.')
                        sign = 0
                    row['sign'] = sign
                rows.append(row)
        df = pd.DataFrame.from_dict(rows)
        return df
