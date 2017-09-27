from .model_checker import stmt_from_rule

def stmts_from_path(path, model, stmts):
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

