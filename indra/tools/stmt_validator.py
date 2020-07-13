class StatementValidator:
    def __init__(self):



class DbRefsEntryValidator:
    @staticmethod
    def validate(entry):
        raise NotImplementedError()


class ChebiPrefix(DbRefsEntryValidator):
    @staticmethod
    def validate(entry):
        return not entry or entry.startswith('CHEBI')


class UniProtIDNotList(DbRefsEntryValidator):
    @staticmethod
    def validate(entry):
        if not isinstance(entry, str):
            return False
        if ',' in entry:
            return False
        return True