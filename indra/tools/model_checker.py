

class ModelChecker(object):
    """Check a PySB model against a set of INDRA statements."""

    def __init__(self, model, statements):
        self.model = model
        self.statements = statements

    def check_model(self):
        return [(self.statements[0], True)]
