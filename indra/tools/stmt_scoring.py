from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import pickle
import random
from indra.statements import *

class StmtScoring(object):

    def __init__(self):
        self.unscored_stmts = []
        self.scored_stmts = []

    def load_unscored_stmts_from_file(self, stmt_pkl_file):
        """Load a list of statements to score from a pickle file."""
        with open(stmt_pkl_file, 'rb') as f:
            stmts = pickle.load(f)
            self.load_unscored_stmts(stmts)

    def load_unscored_stmts(self, stmts):
        self.unscored_stmts = stmts

    def load_scored_stmts_from_file(self, stmt_pkl_file):
        """Load a list of scored (or partially scored) stmts from pickle file."""
        with open(stmt_pkl_file, 'rb') as f:
            stmts = pickle.load(f)
            self.load_scored_stmts(stmts)

    def load_scored_stmts(self, stmts):
        self.scored_stmts = stmts

    def shuffle_unscored_stmts(self):
        random.shuffle(self.unscored_stmts)

    def begin_scoring(self):
        if not self.unscored_stmts and not self.scored_stmts:
            print("Please load a set of unscored or partially scored "
                  "statements first.")
            return
        stmts_to_score = None
        # If we have unscored statements and no partially scored statements
        # Initialize list of shuffled statements along with empty scoring info.
        if self.unscored_stmts and not self.scored_stmts:
            stmts_to_score = self.unscored_stmts
            already_scored_stmts = []
        # If we have a set of partially scored statements, then score only those
        # with no scoring information yet
        elif self.scored_stmts:
            stmts_to_score = [stmt for stmt, score in self.scored_stmts
                              if score is None]
            already_scored_stmts = [(stmt, score)
                                    for (stmt, score) in self.scored_stmts
                                    if score is not None]
        print("Already:", already_scored_stmts)
        # Sanity check
        assert stmts_to_score is not None
        # If there are no statements to score, that might mean that we've
        # scored them all!
        if stmts_to_score == []:
            print("No statements to score--perhaps all have been scored?")
            return

        # Now, begin scoring session
        for stmt in stmts_to_score:
            print("---------------------------------------------------")
            print("Statement: %s\n" % stmt)
            for ev_ix, ev in enumerate(stmt.evidence):
                print("Evidence %s: %s\n" % (ev_ix, ev.text))
            while True:
                score_input = raw_input('Score (0, 1, q)> ')
                if score_input in ('0', '1'):
                    score_input = int(score_input)
                    break
                elif score_input == 'q':
                    break
                else:
                    print("Please enter 0 or 1.")

            if score_input == 'q':
                break
            else:
                already_scored_stmts.append((stmt, score_input))
                print("Already: %s" % already_scored_stmts)
        # Update self.scored_stmts with new scored statements
        self.scored_stmts = already_scored_stmts
        scored_set = [stmt for stmt, score in already_scored_stmts]
        for stmt in stmts_to_score:
            if stmt not in scored_set:
                self.scored_stmts.append((stmt, None))

        print("\nResulting stmts")
        print('\n'.join([str(s) for s in stmts_to_score]))
        print()
        print('\n'.join([str(s) for s in already_scored_stmts]))
        print()
        print('\n'.join([str(s) for s in self.scored_stmts]))

if __name__ == '__main__':
    a = Agent('A')
    b = Agent('B')
    ev_list = [Evidence(text='asdf'), Evidence(text='qwerty')]
    s1 = Phosphorylation(None, b, evidence=ev_list)
    s2 = Phosphorylation(a, b, evidence=ev_list)
    s3 = Phosphorylation(a, b, 'T', '173', evidence=ev_list)
    s4 = Acetylation(a, b, evidence=ev_list)
    scorer = StmtScoring()
    scorer.load_unscored_stmts([s1, s2, s3, s4])
    scorer.shuffle_unscored_stmts()
    scorer.begin_scoring()
