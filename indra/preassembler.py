import itertools
from indra.statements import *

class Preassembler(object):

    def __init__(self, stmts=None):
        if stmts:
            self.stmts = stmts
        else:
            self.stmts = []

    def add_statements(self, stmts):
        """Add to the current list of statements."""
        self.stmts += stmts

    @staticmethod
    def combine_duplicates(stmts):
        """Combine evidence from duplicate Statements.

        Statements are deemed to be duplicates if they have the same key
        returned by the matches_key() method of the Statement class.

        This function keeps the first instance of each set of duplicate
        statements and merges the evidence lists from all of the other
        statements.

        Returns a list of unique Statements.
        """
        unique_stmts = []
        # Group statements according to whether they are matches (differing
        # only in their evidence)
        for key, duplicates in itertools.groupby(stmts,
                                                 key=lambda x: x.matches_key()):
            # Get the first statement and add the evidence of all subsequent
            # Statements to it
            for stmt_ix, stmt in enumerate(duplicates):
                if stmt_ix == 0:
                    first_stmt = stmt
                else:
                    first_stmt.evidence += stmt.evidence
            # This should never be None or anything else
            assert isinstance(first_stmt, Statement)
            unique_stmts.append(first_stmt)
        return unique_stmts

    def assemble(self):
        self.unique_stmts = Preassembler.combine_duplicates(self.stmts)

        """
        for s1 in self.stmts:
            duplicate_list = []
            for s2 in self.stmts:
                if s1 == s2:
                    continue
                elif s1.matches(s2):
                    s1.supported_by.append(s2)
                    s2.supports.append(s1)
                    duplicate_list.append(s2)
            for dup in duplicate_list:
                self.stmts.remove(dup)
        # Next, refinement/subsumption relationships are considered.
        # Note that since a given statement could potentially support multiple
        # parent statements (e.g., RAF phosphorylates MEK), statements can't
        # be thrown away after the first relationship is found.
        print "self.stmts", self.stmts
        # A list of groups of statements. Each entry is a group (a list of
        # statements of the same type, involving the same Agents).
        stmt_groups = []
        for s1, s2 in itertools.combinations(self.stmts, 2):
            # If the two statements are of the same type and deal with the same
            # Agents (as determined by name or grounding), then we check
            # refinement relationships.
            if refines(s1, s2):
                pass

    def get_statement_family(stmt):
        pass
        """


if __name__ == '__main__':
    #import pickle
    #with open('example_kami.pck') as f:
    #    stmts = pickle.load(f)
    raf = Agent('RAF1')
    mek = Agent('MEK1')
    erk = Agent('ERK2')
    p1 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='foo'))
    p2 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='bar'))
    p3 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='baz'))
    p4 = Phosphorylation(raf, mek, 'Phosphorylation', None,
            evidence=Evidence(text='beep'))
    p5 = Phosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='foo'))
    p6 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='bar'))
    p7 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='baz'))
    p8 = Dephosphorylation(mek, erk, 'Phosphorylation', None,
            evidence=Evidence(text='beep'))
    p9 = Dephosphorylation(Agent('SRC'), Agent('KRAS'),
                         'Phosphorylation', None, evidence=Evidence(text='beep'))
    stmts = [p1, p2, p3, p4, p5, p6, p7, p8, p9]
    """
    for group in itertools.groupby(stmts, key=lambda s: type(s).__name__):
        print '----', group[0], '------'
        for thing in group[1]:
            print thing
    """
    pa = Preassembler(stmts)
    pa.assemble()

    """
    for s in stmts:
        print '-----------'
        s.print_supports()
    """
    # Two duplicate statements, both have strings, evidence should be merged
    # into lists
    # Two dupes, source has list, other has string, string is added to list
    # Two dupes source has string, other has list, list is merged
    # Two dupes, both lists, lists are merged
    # All of the above, check citations, evidence, annotations, etc.
    # What if cell type is different???

    # Add supported_by field (can be more than one obj)
    # Add supports field (can be more than one obj)

    # >2 duplicate statements, should end up with one top level, with multiple
    # supported_by (not a stack)?
