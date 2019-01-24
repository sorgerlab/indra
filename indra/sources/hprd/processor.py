from indra.statements import Complex, Agent
from indra.sources.hprd.api import HprdComplexRow

class HprdProcessor(object):
    def __init__(self, hprd_rows):
        hprd_complex_rows = [r for r in hprd_rows
                             if isinstance(r, HprdComplexRow)]
        self.statements = []
        self.get_complexes(hprd_complex_rows)

    def get_complexes(self, hprd_complex_rows):
        stmt = Complex([Agent('ASCL1'), Agent('TCF3'), Agent('MEF2C')])
        self.statements.append(stmt)
