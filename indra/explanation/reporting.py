from indra.sources.indra_db_rest.api import get_statements_by_hash


def stmts_from_pysb_path(path, model, stmts):
    """Return source Statements corresponding to a path in a model.

    Parameters
    ----------
    path : list[tuple[str, int]]
        A list of tuples where the first element of the tuple is the
        name of a rule, and the second is the associated polarity along
        a path.
    model : pysb.core.Model
        A PySB model which contains the rules along the path.
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements from which the model was assembled.

    Returns
    -------
    path_stmts : list[indra.statements.Statement]
        The Statements from which the rules along the path were obtained.
    """
    path_stmts = []
    for path_rule, sign in path:
        for rule in model.rules:
            if rule.name == path_rule:
                stmt = stmt_from_rule(path_rule, model, stmts)
                assert stmt is not None
                path_stmts.append(stmt)
    return path_stmts


def stmts_from_indranet_path(path, model, signed, from_db=True, stmts=None):
    """Return source Statements corresponding to a path in an IndraNet model
    (found by SignedGraphModelChecker or UnsignedGraphModelChecker).

    Parameters
    ----------
    path : list[tuple[str, int]]
        A list of tuples where the first element of the tuple is the
        name of an agent, and the second is the associated polarity along
        a path.
    model : nx.Digraph or nx.MultiDiGraph
        An IndraNet model flattened into an unsigned DiGraph or signed
        MultiDiGraph.
    signed : bool
        Whether the model and path are signed.
    from_db : bool
        If True, uses statement hashes to query the database. Otherwise, looks
        for path statements in provided stmts.
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements from which the model was assembled. Required
        if from_db is set to False.

    Returns
    -------
    path_stmts : list[[indra.statements.Statement]]
        A list of lists of INDRA statements explaining the path (each inner
        corresponds to one step in the path because the flattened model can
        have multiple statements per edge).
    """
    steps = []
    for i in range(len(path[:-1])):
        source = path[i]
        target = path[i+1]
        if signed:
            if source[1] == target[1]:
                sign = 0
            else:
                sign = 1
            stmt_data = model[source[0]][target[0]][sign]['statements']
        else:
            stmt_data = model[source[0]][target[0]]['statements']
        hashes = [stmt['stmt_hash'] for stmt in stmt_data]
        if from_db:
            statements = get_statements_by_hash(hashes, simple_response=True)
        else:
            statements = [
                stmt for stmt in stmts if stmt.get_hash() in hashes]
        steps.append(statements)
    return steps


def stmts_from_pybel_path(path, model, from_db=True, stmts=None):
    """Return source Statements corresponding to a path in a PyBEL model.

    Parameters
    ----------
    path : list[tuple[str, int]]
        A list of tuples where the first element of the tuple is the
        name of an agent, and the second is the associated polarity along
        a path.
    model : pybel.BELGraph
        A PyBEL BELGraph model.
    from_db : bool
        If True, uses statement hashes to query the database. Otherwise, looks
        for path statements in provided stmts.
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements from which the model was assembled. Required
        if from_db is set to False.

    Returns
    -------
    path_stmts : list[[indra.statements.Statement]]
        A list of lists of INDRA statements explaining the path (each inner
        corresponds to one step in the path because PyBEL model can have
        multiple edges representing multiple statements and evidences between
        two nodes).
    """
    steps = []
    for i in range(len(path[:-1])):
        source = path[i]
        target = path[i+1]
        edges = model[source[0]][target[0]]
        hashes = set()
        for j in range(len(edges)):
            try:
                hashes.add(edges[j]['stmt_hash'])
            # If a statement subject or object is a Complex, model would have
            # edges with non-regular indices
            except KeyError:
                continue
        if from_db:
            statements = get_statements_by_hash(list(hashes),
                                                simple_response=True)
        else:
            statements = [
                stmt for stmt in stmts if stmt.get_hash() in hashes]
        steps.append(statements)
    return steps


def stmt_from_rule(rule_name, model, stmts):
    """Return the source INDRA Statement corresponding to a rule in a model.

    Parameters
    ----------
    rule_name : str
        The name of a rule in the given PySB model.
    model : pysb.core.Model
        A PySB model which contains the given rule.
    stmts : list[indra.statements.Statement]
        A list of INDRA Statements from which the model was assembled.

    Returns
    -------
    stmt : indra.statements.Statement
        The Statement from which the given rule in the model was obtained.
    """
    stmt_uuid = None
    for ann in model.annotations:
        if ann.subject == rule_name:
            if ann.predicate == 'from_indra_statement':
                stmt_uuid = ann.object
                break
    if stmt_uuid:
        for stmt in stmts:
            if stmt.uuid == stmt_uuid:
                return stmt
